// To disable the context menu (for HTML inspecting), uncomment the line with the @@ comment.

// Apply edit mode to buttons on load
if (window.location.hash == '#edit') {
	toggleEditModeBtns(true)
}

// Prevent accidentally closing or leaving the page.
window.onbeforeunload = e => {
	if (table.isSubmitting) {
		table.isSubmitting = false
	} else if (table.hasEdits()) {
		return confirm('Are you sure?')
	}
}

/////////////////////////////////
// #region - Table rendering

// Parse data
const dataInput = JSON.parse(document.getElementById('table').getAttribute('data'))
document.getElementById('table').removeAttribute('data')

// Assemble context menus
const contextMenus = new ContextMenus()

// Parse columns
const columns = parseColumns(dataInput)

// Find or create a unique index column
const { index, addedIndex } = ensureUniqueIndex(dataInput, columns)

// Create table
const table = new Table('#table', {
	// Tabulator options
	data: dataInput,
	columns,
	index,
	pagination: false,
	editMode: window.location.hash == '#edit',
	movableRows: false, // We use our own custom row reordering - see table.onCellMouseDown()
	movableColumns: true,
	rowRange: 'active', // Master checkbox will only select filtered rows

	// Non-Tabulator options
	addedIndex,

	// We're not using the built-in keybindings, as
	// we need more advanced logic. See 'keydown' event listener.
	keybindings: false,

	// We're not using reactive data because it's
	// not reflecting deleted rows and thus useless.
	// Instead we run table.getDataFinal().
	// - - -
	// reactiveData: true,

	// We're not using the built-in selection
	// mechanism – see onRowClick()
	// - - -
	// selectable: true,
	// selectableRangeMode: 'click',

	// We're not letting users resize rows,
	// because it complicates truncation, and
	// also because we can't reset the row height
	// the same way as we do with the columns.
	// - - -
	// resizableRows: false,
})

// Pass table & action references to the context menu class.
contextMenus.init(table, {
	saveEdits: () => {
		toggleEditMode(false)
	},
	cancelEdits: () => {
		toggleEditMode(false, true)
	},
})

// Display action dropdown.
table.on('rowSelected', () => {
	toggleSelectionActions(true)
})

// Hide action dropdown
table.on('rowDeselected', () => {
	if (!table.getSelectedRows().length) {
		toggleSelectionActions(false)
	}
})

//
//

// Create the columns object that is fed to Tabulator.
function parseColumns(data) {
	// Get cell properties for each column.
	const cellProps = _getCellProps(data)

	// Get sorters and formatters for each column.
	const { sortersByColumn, formattersByColumn } = _getSortersAndFormatters(cellProps)

	// Create columns
	const columns = Object.entries(cellProps).map(([colName, propsPerRow]) => {
		const { formatter, formatterParams } = formattersByColumn[colName]

		col = {
			sorter: sortersByColumn[colName],
			title: colName,
			field: colName,
			width: propsPerRow.some(cellProps => cellProps.isLongText) ? 400 : undefined, // Set column width only for columns that contain long text,
			editor: _customTextareaEditor,
			editorParams: { selectContents: true }, // Select value on focus,
			editable: () => table.isEditMode(),
			formatter,
			formatterParams,
			headerMenu: contextMenus.header,
			headerContextMenu: contextMenus.header,
			contextMenu: contextMenus.cell, // @@
		}

		// Add sorter params
		if (propsPerRow == 'date') {
			col.sorterParams = {
				format: 'yyyy-MM-dd',
				alignEmptyValues: 'top',
			}
		} else if (propsPerRow == 'number') {
			col.sorterParams = {
				thousandSeparator: ',',
				decimalSeparator: '.',
				alignEmptyValues: 'top',
			}
		}

		return col
	})

	// For debugging
	// - - -
	// console.log('Sorters:', sorters)
	// console.log('Data:', data)
	// console.log('Columns:', columns)
	// console.log('ContentProps', contentPropsAll)
	// console.log('needFormatter', needFormatter)
	// console.log(cellProps)
	// console.log('formattersByColumn:', formattersByColumn)
	return columns
}

// Get sorters and formatters for each column.
function _getSortersAndFormatters(cellProps) {
	const sortersByColumn = {}
	const formattersByColumn = {}

	Object.entries(cellProps).forEach(([colName, propsPerRow]) => {
		// Select sorter
		// - - -
		// To figure out the correct column sorter, we check if
		// the type is the same for all rows (ignore blank values),
		// and if it is, we set it as the column's sorter.
		// - - -
		// Supported sorters can be found here:
		// https://tabulator.info/docs/5.5/sort#func-builtin
		if (propsPerRow.every((cellProps, i, arr) => cellProps.type == arr[0].type || cellProps.isEmpty)) {
			sortersByColumn[colName] = propsPerRow[0].type
		} else {
			sortersByColumn[colName] = 'string'
		}

		// Select formatter
		// - - -
		// To figure out the formatter, we check if all rows in
		// a particular column have a certain property (like isUrl)
		// and use the appropriate formatter if they do (igoring empty
		// values). All the rest uses the default formatter.
		if (propsPerRow.every(cellProps => cellProps.isUrl || cellProps.isEmpty)) {
			// Url formatter (built-in)
			formattersByColumn[colName] = {
				formatter: 'link',
				formatterParams: {
					target: '_blank',
					label: cell => {
						let val = cell.getValue()
						return val ? val.replace(/^http(s)?:\/\/([a-zA-Z0-9$-_.+!*'(),/&?=:%]+?)(\/)?$/, '$2') : val // Reformat URLs
					},
				},
			}
		} else {
			// Text formatter (custom, letting us use advanced truncation)
			formattersByColumn[colName] = {
				// This custom formatter lets us truncate text to a custom number
				// of lines. It's in lieue of the built-un 'textarea' formatter.
				formatter: (cell, formatterParams, onRendered) => {
					return `<div class="text-wrap">${cell.getValue()}</div>`
				},
				formatterParams: null,
			}
		}
	})

	return { sortersByColumn, formattersByColumn }
}

// Returns an object with arrays of cell property objects per row:
// { col1: [
// 	  { type: 'string', isUrl: false, isLongText: false }
// 	  { type: 'string', isUrl: true, isLongText: false }
// 	  { type: 'string', isUrl: false, isLongText: true }
//   ],
//   col2: ...
// }
function _getCellProps(data) {
	const cellProps = {}

	data.forEach((row, i) => {
		Object.entries(row).map(([colName, val]) => {
			val_str = val ? val.toString() : ''

			if (!cellProps[colName]) cellProps[colName] = []
			cellProps[colName].push({})

			// Set content type.
			if (isDate(val_str)) {
				cellProps[colName][i].type = 'date'
				row[colName] = moment(val_str).format('YYYY-MM-DD') // Reformat dates
			} else {
				const type = val === null ? 'string' : typeof val
				cellProps[colName][i].type = typeof val
			}

			// Set isUrl
			cellProps[colName][i].isUrl = Boolean(val_str.match(/^http(s)?:\/\//))
			// row[key] = val.replace(/^http(s)?:\/\/([a-zA-Z0-9$-_.+!*'(),/&?=:%]+?)(\/)?$/, '$2') // Reformat URLs

			// Set isLongText
			cellProps[colName][i].isLongText = val_str.length > 70

			// Set isEmpty
			cellProps[colName][i].isEmpty = val === null
		})
	})
	return cellProps
}

// There's a built-in 'textarea' editor, but we weren't happy with
// how truncation was handled for cells containing long text, so
// we had to create a custom formatter that gives us more control
// (specifically that lets us truncate text to a custom number of
// lines) and a custom formatter begets a custom editor.
// - - -
// About custom editors: https://tabulator.info/docs/5.5/edit#edit-custom
// About built-in editors: https://tabulator.info/docs/5.5/edit#edit-builtin
function _customTextareaEditor(cell, onRendered, success, cancel, editorParams) {
	const editor = document.createElement('textarea')
	editor.value = cell.getValue()

	// This ensures the textarea height will not exceed the number of lines.
	editor.style.overflow = 'hidden'

	onRendered(() => {
		// Expand textarea to fit content
		editor.style.height = ''
		editor.style.height = editor.scrollHeight + 4 + 'px' // +4px to make up for border
		editor.style.overflow = 'auto'

		// Focus or select text
		if (editorParams.selectContents) {
			editor.select()
		} else {
			editor.focus()
		}
	})

	// Note: blur and change usually fire together, but we can't
	// limit it to just change, becayse
	// editor.addEventListener('change', _applyEdit)
	editor.addEventListener('blur', _applyEdit)
	editor.addEventListener('keydown', _onKeyDown)
	editor.addEventListener('input', _onInput)

	return editor

	//
	//

	function _onKeyDown(e) {
		if (e.key == 'Enter') {
			// This is the same behavior as a Google sheet:
			// Enter will exit the edit field, but cmd + enter
			// will result in line break. But because with the meta
			// key pressed, there is no actual character written,
			// we have to add the line break ourselves.
			if (e.metaKey || e.ctrlKey) {
				editor.value += '\n'
				_onInput()
			} else {
				_applyEdit()
			}
		} else if (e.key == 'Escape') {
			cancel()
		}
	}

	// Resize the textarea after the character was registered.
	function _onInput() {
		editor.style.height = editor.scrollHeight + 4 + 'px'
	}

	function _applyEdit() {
		success(editor.value)
	}
}

// #endregion

/////////////////////////////////
// #region - Data prepping

// Make sure our data has a unique index.
// - - -
// When updating content on save, we depend on updateData()
// because setData() causes the entire table to redraw and
// the page to scroll up. But updateData relies on an index
// to be defined, i.e. a column with unique values. Because
// we don't know what kind of data will be loaded, we need to
// do some trickery. In this order:
// 1. Check if there already is a unique index column
// 2. If there is, check if the values are unique
// 3. If no valid index column is found, we create our own '#'
function ensureUniqueIndex(data, columns) {
	let index = _findUniqueIndex(data, columns)
	let addedIndex = false
	const $dropdownExportIndex = document.getElementById('opt-export-index-col')
	if (!index) {
		// Add index column
		const indexName = _addIndexCol(data, columns)
		index = indexName
		addedIndex = true
	} else {
		// Use existing index column

		// Set the default value for the "export index column" option in the options panel.
		// Note: this runs before carbon.js has loaded, otherwise we'd have to update the
		// value using carbonUi.updateDropdown('#opt-export-index-col', '1')
		$dropdownExportIndex.value = '1'

		// Make index column non-editable,remove menu and sent to the front.
		const indexCol = columns.find(col => col.field == index)
		indexCol.editable = false
		indexCol.headerMenu = false
		columns.splice(columns.indexOf(indexCol), 1)
		columns.unshift(indexCol)
	}

	// Update dropdown label to reflect the index column's name.
	const $dropdownLabel = $dropdownExportIndex.parentElement.querySelector('.ibm-label-txt')
	$dropdownLabel.innerHTML = $dropdownLabel.innerHTML.replace(/index/, `'${index}'`)

	return { index, addedIndex }
}

// Check if there already is a unique index column
function _findUniqueIndex(data, columns) {
	const fields = columns.map(col => col['field'].toLowerCase())

	// Check for 'index' column (case-insensitive)
	if (fields.includes('index')) {
		const field = columns[fields.indexOf('index')]['field']
		if (_areColumnValuesIncrementalIndex(field, data)) {
			return field
		}
	}

	// Check for '#' column
	if (fields.includes('#')) {
		if (_areColumnValuesIncrementalIndex('#', data)) {
			return '#'
		}
	}

	return false
}

// Check if the values of a suspected index column are incremental.
// Can start from any number but needs to be consecutive to be valid.
function _areColumnValuesIncrementalIndex(field, data) {
	// return false // For testing $$
	const values = data.map(row => +row[field])
	const valuesSorted = [...values].sort((a, b) => a - b)
	const firstValue = +valuesSorted[0]
	if (typeof firstValue != 'number') return false
	return valuesSorted.every((val, i) => val == i + firstValue)
}

// Add index column if none is found.
// This will usually be '#' but we gotta make sure there's no conflict.
function _addIndexCol(data, columns) {
	// Find a unique name for the index column that is not already taken.
	const fields = columns.map(col => col['field'].toLowerCase())
	const indexNameOptions = ['#', 'index', 'idx', 'nr', '*', '-']
	let i = 0
	let indexName = indexNameOptions[i]
	while (fields.includes(indexName)) {
		i++
		indexName = indexNameOptions[i]
	}

	// Add index column.
	data.forEach((row, i) => {
		row[indexName] = i + 1
	})

	columns.unshift({
		sorter: 'number',
		title: indexName,
		field: indexName,
		headerContextMenu: contextMenus.header,
		contextMenu: contextMenus.cell,
	})

	return indexName
}

//#endregion

/////////////////////////////////
// #region - Event listeners

// Edit
document.getElementById('btn-edit').addEventListener('click', () => {
	toggleEditMode(true)
})

// Cancel
document.querySelector('#btn-cancel').addEventListener('click', () => {
	toggleEditMode(false, true)
})

// Done
document.querySelector('#btn-done').addEventListener('click', () => {
	toggleEditMode(false)
})

// Submit
document.querySelector('#btn-submit').addEventListener('click', submitData)

// Deselect
document.querySelector('#btn-deselect').addEventListener('click', table.deselectRows.bind(table))

// Selection actions
document.querySelector('#dd-selection-actions').addEventListener('change', dispatchSelectionActions)

// General actions
document.querySelector('#dd-actions').addEventListener('change', dispatchMainActions)

// Display options
document.querySelector('#btn-options').addEventListener('click', () => {
	toggleOptions(true)
})

// Submit options
document.querySelector('#options .btn-apply').addEventListener('click', () => {
	applyOptions()
	toggleOptions(false)
})

// Cancel options
document.querySelector('#options .btn-cancel').addEventListener('click', () => {
	toggleOptions(false)
})

// Reset columns
document.querySelector('#reset-links .reset-col-width').addEventListener('click', e => {
	table.resetColWidths()
	e.preventDefault()
})

// Header click --> scroll to top
document.getElementById('cloak').addEventListener('click', () => {
	window.scrollTo(0, 0)
})

// Exit focus on blur
document.addEventListener('click', e => {
	if (e.target.closest('.tabulator-cell') || e.target.classList.contains('tabulator-cell')) return
	table.unsetFocus()
})

// #endregion

/////////////////////////////////
// #region - Option panel

// Toggle options panel
function toggleOptions(bool) {
	bool = bool != undefined ? bool : !document.getElementById('options').classList.contains('show')
	if (bool) {
		document.getElementById('options').classList.remove('hide')
		document.getElementById('btn-options').classList.add('hide')
	} else {
		document.getElementById('options').classList.add('hide')
		document.getElementById('btn-options').classList.remove('hide')
	}
}

// Apply options
function applyOptions() {
	// Truncation
	const truncLimit = +document.getElementById('dd-opt-select-truncation').value
	_setTruncation(truncLimit)
}

// Set truncation limit (0/1/3)
function _setTruncation(limit) {
	const $table = document.getElementById('table')
	$table.classList.remove('trunc-3', 'trunc-1')
	if (limit > 0) {
		$table.classList.add('trunc-' + limit)
	}
	table.redraw(true)
}

// #endregion

/////////////////////////////////
// #region - Utility functions

// Toggle table edit mode
function toggleEditMode(bool, revertChanges) {
	table.toggleEditMode(bool, revertChanges)
	toggleEditModeBtns(bool)
}

// Toggle buttons edit mode
function toggleEditModeBtns(bool) {
	if (bool) {
		document.getElementById('btn-wrap-left').classList.add('edit-mode')
	} else {
		document.getElementById('btn-wrap-left').classList.remove('edit-mode')
	}
}

// Toggle the dropdown with actions relating to your current selection.
function toggleSelectionActions(bool) {
	if (bool) {
		document.getElementById('btn-wrap-left').classList.add('show-actions')
	} else {
		document.getElementById('btn-wrap-left').classList.remove('show-actions')
	}
}

// Check if string is a date.
function isDate(str, log) {
	// Ignore numbers that are not formatted like a date.
	const separators = str.match(/[-./]/g, '')
	if (separators && separators.length == 2 && separators[0] == separators[1]) {
		return moment(str).isValid()
	} else if (str.match(/(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)/i)) {
		return moment(str).isValid()
	} else {
		return false
	}
}

// Dispatch actions from the selection dropdown.
function dispatchSelectionActions(e) {
	const action = e.target.value
	if (action == 'delete') {
		contextMenus.deleteSelected()
	} else if (action == 'keep') {
		contextMenus.keepSelected()
	} else if (action == 'copy') {
		contextMenus.copyData(true)
	} else if (action == 'download') {
		contextMenus.downloadData(true)
	}
	e.target.value = ''
}

// Dispatch actions from the action dropdown.
function dispatchMainActions(e) {
	const action = e.target.value
	if (action == 'copy') {
		contextMenus.copyData()
	} else if (action == 'download') {
		contextMenus.downloadData()

		// Note: we're not using the built-in download function because
		// we built our own movableRows UI, see table.onCellMouseDown().
		// table.download('csv', 'data.csv')
	}
	e.target.value = ''
}

// #endregion

/////////////////////////////////
// #region - Submit data

// Submit data back to CLI/Jupyter
function submitData() {
	table.isSubmitting = true
	const data = table.getDataFinal()

	// Create a new XMLHttpRequest object
	var xhr = new XMLHttpRequest()

	// Define the request method and URL.
	xhr.open('POST', '/submit', true)

	// Set up a callback function to handle the response.
	xhr.onload = function () {
		if (xhr.status === 200) {
			// Success
			window.location.href = `/success` // ?data=${data}
		} else {
			// Error
			alert('Submit request failed with status code ' + xhr.status)
		}
	}

	// Send the request
	xhr.send(JSON.stringify(data))
}

//#endregion
