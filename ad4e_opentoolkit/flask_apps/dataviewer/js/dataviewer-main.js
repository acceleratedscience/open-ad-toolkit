// To disable the context menu (for HTML inspecting), uncomment the line with the @@ comment.

// Apply edit mode to buttons on load
if (window.location.hash == '#edit') {
	toggleEditModeBtns(true)
}

// Prevent accidentally closing or leaving the page.
window.onbeforeunload = e => {
	if (table.hasEdits()) return confirm('Are you sure?')
}

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
document.querySelector('#btn-deselect').addEventListener('click', deselectRows)

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
	table.resetCols()
	e.preventDefault()
})

// Header click --> scroll to top
document.getElementById('cloak').addEventListener('click', () => {
	window.scrollTo(0, 0)
})

// document.querySelector('#reset-links .reset-row-width').addEventListener('click', e => {
// 	table.resetRows()
// 	e.preventDefault()
// })

// document.addEventListener('keydown', e => {
// 	e.preventDefault()
// })

// Key handlers
document.addEventListener('keydown', e => {
	const inputInFocus = e.target.tagName.toLowerCase() == 'input' || e.target.tagName.toLowerCase() == 'textarea'
	if (inputInFocus) return

	const $focusCell = document.querySelector('.tabulator-cell.focus')
	const $row = $focusCell ? $focusCell.closest('.tabulator-row') : null
	const row = $row ? table.getRow($row) : null

	// Copy focused cell content on cmd + C
	if ($focusCell && e.key == 'c' && (e.metaKey || e.ctrlKey)) {
		contextMenus.copyCell(null, null, $focusCell)
		return
	}

	if (e.key == 'Escape') {
		// Hide full cell content overlay
		const fullCellDisplayActive = hideFullCellContent()
		if (fullCellDisplayActive) return

		// Deselect rows
		if (table.getSelectedRows().length) {
			deselectRows(true)
			return
		}

		// Exit edit mode
		if (table.isEditMode()) {
			toggleEditMode(false, true)
			return
		}
	} else if ((e.metaKey || e.ctrlKey) && e.key == 'z') {
		// Undo / redo
		if (e.shiftKey) {
			table.redo()
		} else {
			table.undo()
		}
		return
	}

	// Move focus with arrow keys
	if ($focusCell) {
		if (e.key == 'ArrowLeft' || (e.shiftKey && e.key == 'Tab')) {
			moveFocus('left')
			e.preventDefault()
		} else if (e.key == 'ArrowRight' || e.key == 'Tab') {
			moveFocus('right')
			e.preventDefault()
		} else if (e.key == 'ArrowUp') {
			moveFocus('up')
			e.preventDefault()
		} else if (e.key == 'ArrowDown') {
			moveFocus('down')
			e.preventDefault()
		} else if (e.key == 'Escape') {
			$focusCell.classList.remove('focus')
			e.preventDefault()
		} else if (e.key == 'Enter') {
			if (table.isEditMode()) {
				// Edit mode --> edit
				setTimeout(() => {
					$focusCell.click()
				}, 0)
			} else {
				// View mode
				if (e.metaKey || e.ctrlKey) {
					// meta --> edit
					toggleEditMode(true)
					setTimeout(() => {
						$focusCell.click()
					}, 0)
				} else {
					// --> select
					row.toggleSelect()
				}
			}
		} else if (table.isEditMode() && e.key.length == 1) {
			// Edit mode --> start typing
			toggleEditMode(true)
			setTimeout(() => {
				$focusCell.click()
				const $input = $focusCell.querySelector('input') || $focusCell.querySelector('textarea')
				if ($input) {
					$input.value = e.key
				}
			}, 0)
		}
	}
})

// Exit focus on blur
document.addEventListener('click', e => {
	if (e.target.closest('.tabulator-cell') || e.target.classList.contains('tabulator-cell')) return
	unsetFocus()
})

// #endregion

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

// Create table ##
const table = new Table('#table', {
	// Tabulator options
	data: dataInput,
	columns,
	index,
	pagination: false,
	editMode: window.location.hash == '#edit',
	movableRows: false,
	movableColumns: true,
	rowRange: 'active', // Master checkbox will only select filtered rows
	rowFormatter: row => {
		// const $row = row.getElement()
		// const rowHeight = $row.clientHeight
		// const lines = Math.round((rowHeight - 8) / 20) // Math round not necessary but just in case
		// // console.log(lines)
		// if (lines == 1) {
		// 	$row.classList.add('trunc-1')
		// }
		// console.log('FORMAT')
	},

	// %% trash
	// rowFormatter: row => {
	// 	const cells = row.getCells()
	// 	cells.forEach(cell => {
	// 		const $cell = cell.getElement()
	// 		const $textWrap = $cell.querySelector('.text-wrap')
	// 		if ($textWrap) {
	// 			$cell.classList.add('long-text')
	// 		}
	// 	})
	// },

	// Disable built-in keybindings, we need more
	// advanced logic. See 'keydown' event listener.
	keybindings: false,

	// Non-Tabulator options
	addedIndex,

	// We're not using reactive data because it's
	// not reflecting deleted rows and thus useless.
	// Instead we run prepDataForExport().
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
	// resizableRows: true,
})

// Pass table reference to context menus
contextMenus.init(table, {
	saveEdits: () => {
		toggleEditMode(false)
	},
	cancelEdits: () => {
		toggleEditMode(false, true)
	},
	deselectRows,
})

table.on('tableBuilt', () => {
	// Redraw the table after the data object is built but before the table is rendered.
	// This is a little hack to prevent all rows to take on the height of the tallest cell.
	table.redraw(true)

	// Store the initial state of the table.
	table.add_history(table.getData())

	// %% trash
	// // Set the "Add index column" dropdown to the correct value.
	// setIndexColDropdown()
})

// Controls our homebrewed movableRows.
table.on('cellMouseDown', (e, row) => {
	table.onCellMouseDown(e, row, onCellClick)
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

// // Set (artificial) focus on edited cell.
// // cellEditing
// table.on('cellEdited', cell => {
// 	$focusCell = cell.getElement()
// 	setFocus($focusCell)
// 	// setTimeout(() => { // %% trash
// 	// 	console.log(22, $focusCell)
// 	// 	setFocus($focusCell)
// 	// }, 1000)
// })

// To do: implement short history so you can undo changes.
table.on('dataChanged', data => {
	table.add_history(data)
	console.log('dataChanged', data)
})

// Update truncation whenever text is changed.
table.on('cellEdited', cell => {
	const index = cell.getRow().getIndex()
	const field = cell.getField()
	const row = table.getRows()[index - 1]
	const $row = row.getElement()
	const $cell = row.getCell(field).getElement()

	// // First remove possible 1-line truncation from the row.
	// $row.classList.remove('trunc-1')

	// This will trigger the rowFormatter, which will
	// re-assign the trun-1 clas unless any of the cells
	// in this row have more than one line of text.
	table.redraw(true)

	// Trace back the same cell after redraw has replaced
	// the HTML and reset de focus.
	const $newCell = table.getRows()[index - 1].getCell(field).getElement()
	setFocus($newCell)
})

// %% trash
// table.on('dataSorted', function (sorters, rows) {
// 	// Ignore when triggered during initiation.
// 	// Can be recoignized by the empty sorters array.
// 	if (sorters.length) table.sortChange = true
// })

// Note: onRowClick is controlled from within onCellClick
// this way we can control how they play together.
// table.on('cellClick', onCellClick)

// // Double click doesn't play well with selection
// table.on('cellDblClick', e => {
// 	$cell = e.target.classList.contains('tabulator-cell') ? e.target : e.target.closest('.tabulator-cell')
// 	copyCellContent($cell)
// })

//
//

// Create the columns object that is fed to Tabulator.
function parseColumns(data) {
	// Get cell properties for each column.
	const cellProps = getCellProps(data)

	// Get sorters and formatters for each column.
	const { sortersByColumn, formattersByColumn } = getSortersAndFormatters(cellProps)

	// console.log(cellProps)
	// console.log('Sorters:', sorters)
	// console.log('formattersByColumn:', formattersByColumn)

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
			// headerSort: false, %% trash
			// headerClick: onHeaderClick, // We use our own sort logic via headerClick - %% trash

			// // This blocks the user from resizing the
			// // column, so we use a formatter instead:
			// maxWidth: 500,
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
	// console.log('Sorters:', sorters)
	// console.log('Data:', data)
	// console.log('Columns:', columns)
	// console.log('ContentProps', contentPropsAll)
	// console.log('needFormatter', needFormatter)
	return columns
}

// Get sorters and formatters for each column.
function getSortersAndFormatters(cellProps) {
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
function getCellProps(data) {
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

// Look at a column's content and pick the appropriate formatter.
// The formatter defines how a cell is rendered.
// - - -
// Returns { formatter, formatterParams, formatterType }
function _pickFormatter(key, needFormatter, contentPropsAll) {
	const hasSomeLongText = contentPropsAll[key].some(props => props.isLongText)
	const isAllUrls = contentPropsAll[key].every(props => props.isUrl)

	if (needFormatter[key]) {
		if (isAllUrls) {
			return { formatter: 'link', formatterParams: { target: '_blank' } }
		} else if (hasSomeLongText) {
			return {
				formatterType: 'textarea',
				// This custom formatter lets us truncate text to a custom number
				// of lines. It's in lieue of the built-un 'textarea' formatter.
				formatter: (cell, formatterParams, onRendered) => {
					return `<div class="text-wrap">${cell.getValue()}</div>`
				},
			}
		}
	}
	return {}
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

	onRendered(() => {
		// Expand textarea to fit content
		editor.style.height = ''
		editor.style.height = editor.scrollHeight + 4 + 'px' // +4px to make up for border

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
		console.log('#')
		editor.style.height = editor.scrollHeight + 4 + 'px'
	}

	function _applyEdit() {
		console.log('_applyEdit')
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

// %% trash
// // Check if the values of a certain column are unique
// function _areColumnValuesUnique(field, data) {
// 	for (let i in data) {
// 		const val = data[i][field]
// 		for (let j in data) {
// 			if (i != j && val == data[j][field]) {
// 				return false
// 			}
// 		}
// 	}
// 	return true
// }

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
// #region - Data modification

// Delete column and all its data
function deleteColumn(e, col) {
	const field = col.getField()
	if (field == table.options.index) {
		let message = `You can't delete the '${field}' column as it's used for identifying the rows.` // Note: repeat
		if (table.addedIndex) {
			message += ' It will be removed when you submit the data, unless you change this in the options panel.'
		} else {
			message += ' But you can choose not to export it in the options panel.'
		}
		alert(message)
	} else {
		col.delete()
	}
}

//#endregion

/////////////////////////////////
// #region - Table interaction

// %% trash
// function onHeaderClick(e, col) {
// 	table.sort(col)
// }

// %% Trash
// table.on('dataSorting', onDataSorting)
// function onDataSorting(sorters) {
// 	console.log('SORT')
// 	return false
// }

function onCellClick(e, cell) {
	// Activate edit mode when cmd-clicking a cell.
	if (e.metaKey || e.ctrlKey) {
		toggleEditMode(true)
		setTimeout(() => {
			$focusCell.click()
		}, 0)
		return
	}

	// Display overlay with full content when cell is truncated.
	if (e.target.classList.contains('text-wrap')) {
		const $textWrap = e.target
		const $cell = $textWrap.closest('.tabulator-cell')
		const isTruncated = _isTruncated($textWrap)
		if (isTruncated) {
			displayFullCellContent($textWrap, $cell)
			// _expandCell($textWrap, $cell, cell)
			return
		}
	}

	// Set focus on clicked cell.
	$focusCell = cell.getElement()
	setFocus($focusCell)

	// Select row
	onRowClick(e, cell.getRow())

	//
	//

	// Unused but may come in handy later:
	// - - -
	// Expand the cell so all content is visible
	function _expandCell($textWrap, $cell, cell) {
		const $row = $textWrap.closest('.tabulator-row')
		const row = cell.getRow()._row

		$textWrap.classList.add('expand')
		$cell.style.removeProperty('height')

		// Set height of all cells in row to match
		$row.childNodes.forEach($siblingCell => {
			if ($siblingCell.classList.contains('tabulator-cell')) {
				$siblingCell.style.setProperty('height', $cell.offsetHeight + 'px')
			}
		})

		// Update the row height in the table object
		row.height = $cell.offsetHeight
		row.heightStyled = $cell.offsetHeight + 'px'
		row.outerHeight = $cell.offsetHeight + 1
		row.manualHeight = true
	}

	// Check whether a cell has truncated text
	function _isTruncated($textWrap) {
		// Create clone
		const $clone = $textWrap.cloneNode()
		$clone.innerHTML = $textWrap.innerHTML
		_copyStyle($textWrap, $clone)

		// Render clone hidden on DOM
		$clone.style.setProperty('position', 'absolute')
		$clone.style.setProperty('top', 0)
		$clone.style.setProperty('left', 0)
		$clone.style.setProperty('visibility', 'hidden')
		document.body.append($clone)

		// Check if clone is wider than original
		const isTruncated = $clone.offsetHeight > $textWrap.offsetHeight

		// Delete clone
		$clone.remove()
		delete $clone

		return isTruncated
	}

	// Copy styles from one element to another
	function _copyStyle($sourceElm, $targetElm) {
		keys = ['width', 'font-size', 'line-height', 'padding']
		const computedStyle = window.getComputedStyle($sourceElm)
		keys.forEach(key => {
			$targetElm.style.setProperty(key, computedStyle.getPropertyValue(key), computedStyle.getPropertyPriority(key))
		})
	}
}

// The Tabulator implementation of cell selection is rather
// sloppy so we're bypassing it with our own implementation.
// More specifically: the default behavior lets you toggle
// cell selection, but you can't shift-click to select multiple
// cells. There's an option { selectableRangeMode: 'click' }
// which provides the multi-select behabvior, but for some
// reason regular toggling is disabled and you need to cmd-click
// a cell to deselect it, and on top of that it's also not
// supporting the selection of multiple groups of cells.
// It may make sense to contribute a fix to the library at some
// point, but we don't have teh time for that now.
function onRowClick(e, row) {
	const currentRowIndex = row.getPosition()
	const lastSelectedRowIndex = table.lastSelectedRowIndex
	if (e.shiftKey && table.lastSelectedRowSelState != null) {
		// Select or deselect all rows between the last selected row and the current row
		const selectedRows = table.getSelectedRows()
		if (selectedRows.length) {
			let lowIndex = Math.min(lastSelectedRowIndex, currentRowIndex)
			let highIndex = Math.max(lastSelectedRowIndex, currentRowIndex)

			// When you select from bottom to top, we gotta include the highIndex
			// When you select from top to bottom, we gotta include the lowIndex
			if (lowIndex != lastSelectedRowIndex) {
				lowIndex -= 1
				highIndex -= 1
			}

			const toBeSelected = table.getRows().slice(lowIndex, highIndex)
			if (table.lastSelectedRowSelState) {
				table.selectRow(toBeSelected)
				table.lastSelectedRowSelState = true
			} else {
				table.deselectRow(toBeSelected)
				table.lastSelectedRowSelState = false
			}
		}
	} else {
		// Toggle single row
		row.toggleSelect()
		table.lastSelectedRowSelState = table.getSelectedRows().includes(row)
	}
	table.lastSelectedRowIndex = currentRowIndex
	table.selectMode = table.getSelectedRows().length > 0
}

// Display overlay div that matches the cell's content and position,
// but that expands to the bottom to display the full, untruncated content.
function displayFullCellContent($textWrap, $cell) {
	const $display = document.createElement('div')
	$display.setAttribute('id', 'display-full-text')
	$display.setAttribute('tabindex', 0)
	$display.innerHTML = $textWrap.innerHTML
	$cell.append($display)
	$display.focus()
	$display.addEventListener('blur', e => {
		// Don't remove display when display itself is clicked.
		if (e.relatedTarget && e.relatedTarget.querySelector('#display-full-text')) {
			e.target.focus()
		} else {
			$display.remove()
		}
	})
}

// Hide the full cell content overlay
function hideFullCellContent() {
	const $display = document.getElementById('display-full-text')
	if ($display) {
		$display.blur()
		return true
	}
	return false
}

// Deselect all rows, but ask user if it's more than 3
function deselectRows(soft) {
	const selectedRows = table.getSelectedRows()
	// Note: == true is on purpose, because we want to ignore the click event
	if (soft == true && selectedRows.length > 3) {
		if (confirm('Are you sure you want to deselect all rows?')) {
			table.deselectRow()
		}
	} else {
		table.deselectRow()
	}
}

// Set artificial cell focus
function setFocus($focusCell) {
	const $currentFocusCell = document.querySelector('.tabulator-cell.focus')
	if ($currentFocusCell) $currentFocusCell.classList.remove('focus')
	$focusCell.classList.add('focus')
}

// Unset artificial cell focus
function unsetFocus() {
	const $focusCell = document.querySelector('.tabulator-cell.focus')
	if ($focusCell) $focusCell.classList.remove('focus')
}

// Move artificial cell focus
function moveFocus(dir) {
	const $currentFocusCell = document.querySelector('.tabulator-cell.focus')
	if (!$currentFocusCell) return
	const index = Array.from($currentFocusCell.parentNode.querySelectorAll('.tabulator-cell')).indexOf($currentFocusCell)
	let $nextFocusCell = null
	if (dir == 'left') {
		$nextFocusCell = $currentFocusCell.parentNode.querySelectorAll('.tabulator-cell')[index - 1]
	} else if (dir == 'right') {
		$nextFocusCell = $currentFocusCell.parentNode.querySelectorAll('.tabulator-cell')[index + 1]
	} else if (dir == 'up') {
		const $prevRow = $currentFocusCell.closest('.tabulator-row').previousElementSibling
		if ($prevRow) {
			$nextFocusCell = $prevRow.querySelectorAll('.tabulator-cell')[index]
		}
	} else if (dir == 'down') {
		const $nextRow = $currentFocusCell.closest('.tabulator-row').nextElementSibling
		if ($nextRow) {
			$nextFocusCell = $nextRow.querySelectorAll('.tabulator-cell')[index]
		}
	}
	if ($nextFocusCell) {
		$currentFocusCell.classList.remove('focus')
		$nextFocusCell.classList.add('focus')
	}
}

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

	// %% trash
	// // Index column
	// const includeIndexCol = +document.getElementById('opt-add-index-col').value
	// if (includeIndexCol) {
	// 	table.addIndexCol(true)
	// } else {
	// 	table.removeIndexCol()
	// }
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
	} else if (str.match(/[jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec]/i)) {
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
		contextMenus.copyData()
	} else if (action == 'download') {
		contextMenus.downloadData()
	}
	e.target.value = ''
}

// Dispatch actions from the action dropdown.
function dispatchMainActions(e) {
	const action = e.target.value
	if (action == 'copy') {
		contextMenus.copyData(true)
	} else if (action == 'download') {
		contextMenus.downloadData(true)
	}
	e.target.value = ''
}

// %% Trash
// // Set the "Add index column" dropdown to the correct value.
// function setIndexColDropdown() {
// 	if (table.addedIndex) {
// 		document.getElementById('opt-add-index-col').querySelector('option[value="1"]').selected = true
// 	} else {
// 		document.getElementById('opt-add-index-col').querySelector('option[value="0"]').selected = true
// 	}
// }

// #endregion

/////////////////////////////////
// #region - Submit data

function submitData() {
	const data = prepDataForExport()

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

// When deleting a column, the data of this column is not deleted.
// So before we submit the data, we have to first filter it by the
// columns that are still present in the table.
function prepDataForExport() {
	let data = table.getData()
	const fields = table.getColumns().map(col => col.getField())

	// Remove deleted column data.
	data.forEach(row => {
		Object.keys(row).forEach(key => {
			if (!fields.includes(key)) {
				delete row[key]
			}
		})
	})

	// Remove index column per the opt-export-index-col options.
	// Note: The default value will be false if the index column
	// was added by us, and true if it was already present.
	const includeIndex = Boolean(+document.getElementById('opt-export-index-col').value)
	if (!includeIndex) {
		const indexName = table.options.index
		data.forEach(row => {
			delete row[indexName]
		})
	}

	// Reorder data to match the order of the columns.
	data = data.map(row => {
		const newRow = {}
		fields.forEach(field => {
			newRow[field] = row[field]
		})
		return newRow
	})

	return data
}

//#endregion
