/////////////////////////////////
// #region - Event listeners

// Apply edit mode to buttons on load
if (window.location.hash == '#edit') {
	toggleEditModeBtns(true)
}

// Edit
document.getElementById('btn-edit').addEventListener('click', () => {
	toggleEditMode(true)
})

// Cancel
document.querySelector('#btn-cancel').addEventListener('click', () => {
	toggleEditMode(false, true)
})

// Save
document.querySelector('#btn-save').addEventListener('click', () => {
	toggleEditMode(false)
})

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
document.querySelector('#table-wrap > a').addEventListener('click', () => {
	table.resetCols()
})

// Key handlers
document.addEventListener('keydown', e => {
	$cell = document.querySelector('.tabulator-cell.focus')
	$row = $cell ? $cell.closest('.tabulator-row') : null
	row = $row ? table.getRow($row) : null

	// const rowIndex = $row ? Array.from($row.parentNode.querySelectorAll('.tabulator-row')).indexOf($row) : null
	// cellName = $cell ? $cell.getAttribute('tabulator-field') : null
	// row = $cell ? $cell.getRow() : null

	// get indef of $row in its parent div 'tabulator-table'

	// Copy focused cell content on cmd + C
	if ($cell && e.key == 'c' && (e.metaKey || e.ctrlKey)) {
		copyCellContent($cell)
		return
	}

	if (e.key == 'Escape') {
		// Hide full cell content overlay
		const fullCellDisplayActive = hideFullCellContent()
		if (fullCellDisplayActive) return

		// Deselect rows
		deselectRows()
	}

	// Move focus with arrow keys
	if (document.querySelector('.tabulator-cell.focus') && !table.isEditMode()) {
		if (e.key == 'ArrowLeft') {
			moveFocus('left')
			e.preventDefault()
		} else if (e.key == 'ArrowRight') {
			moveFocus('right')
			e.preventDefault()
		} else if (e.key == 'ArrowUp') {
			moveFocus('up')
			e.preventDefault()
		} else if (e.key == 'ArrowDown') {
			moveFocus('down')
			e.preventDefault()
		} else if (e.key == 'Escape') {
			document.querySelector('.tabulator-cell.focus').classList.remove('focus')
			e.preventDefault()
		} else if (e.key == 'Enter') {
			console.log(66, row)
			row.toggleSelect()
			// table.selectRow($row)
			// table.toggleEditMode(true)
			// e.preventDefault()
			// setTimeout(() => {
			// 	document.querySelector('.tabulator-cell.focus').click()
			// }, 100)
		}
	}
})

// Exit focus on blur
document.addEventListener('click', e => {
	if (e.target.closest('.tabulator-cell') || e.target.classList.contains('tabulator-cell')) return
	const $currentFocusCell = document.querySelector('.tabulator-cell.focus')
	if ($currentFocusCell) $currentFocusCell.classList.remove('focus')
})

// #endregion

/////////////////////////////////
// #region - Table rendering

// Create context menus
const contextMenus = {
	header: [
		{
			label: 'Delete column',
			action: (e, col) => {
				col.delete()
			},
		},
	],
	row: [
		{
			label: 'Delete row',
			action: (e, row) => {
				row.delete()
			},
		},
		{
			label: 'Select row',
			action: (e, row) => {
				row.select()
			},
		},
	],
}

// Parse data
const data = JSON.parse(document.getElementById('table').getAttribute('data'))
document.getElementById('table').removeAttribute('data')

// Parse columns
const columns = parseColumns(data)

// Add checkbox column for selection UI
columns.unshift({ formatter: 'rowSelection', titleFormatter: 'rowSelection', align: 'center', headerSort: false })

// Create table
const table = new Table('#table', {
	data,
	columns,
	resizableRows: true,
	// rowHeight: 40,
	reactiveData: true,
	pagination: false,
	editMode: window.location.hash == '#edit',
	movableRows: true,
	movableColumns: true,
	// selectable: true,
	// selectableRangeMode: 'click', // We're noy using this – see onRowClick
	rowRange: 'active', // Master checkbox will only select filtered rows

	// paginationSize: 5,
	// rowContextMenu: contextMenus.row, // @@
})

table.on('tableBuilt', () => {
	// Redraw the table after the data object is built but before the table is rendered.
	// This is a little hack to prevent all rows to take on the height of the tallest cell.
	table.redraw(true)

	// Set the "Add index column" dropdown to the correct value.
	setIndexColDropdown()

	// Init right click menu on table header
	// initRightClickMenu()
})

table.on('cellClick', onCellClick)
table.on('rowClick', onRowClick)
// // Double click doesn't play well with selection
// table.on('cellDblClick', e => {
// 	$cell = e.target.classList.contains('tabulator-cell') ? e.target : e.target.closest('.tabulator-cell')
// 	copyCellContent($cell)
// })

// %% trash
// table.on('headerContext', (e, col) => {
// 	e.preventDefault()
// 	showContextMenu(col)
// })

//
//

// Create the columns object that is fed to Tabulator.
function parseColumns(data) {
	// For each column, we store per column:
	// - An array of sorters by row (guessed by content type). --> Defines how column is being sorted (tabular)
	// - An array of content properties by row. --> Defines the preferred editing UI (custom)
	const sortersAll = {
		// col1: ['string', 'string', 'string']
	}
	const contentPropsAll = {
		// col1: [{ isUrl: false, isLongText: false }, { isUrl: false, isLongText: false }],
	}
	data.forEach((row, i) => {
		Object.entries(row).map(([key, val]) => {
			val_str = val ? val.toString() : ''

			// Set sorters
			if (!sortersAll[key]) sortersAll[key] = []
			if (isDate(val_str)) {
				sortersAll[key].push('date')
				row[key] = moment(val_str).format('YYYY-MM-DD') // Reformat dates
			} else {
				sortersAll[key].push(typeof val)
			}

			// Set content props
			if (!contentPropsAll[key]) contentPropsAll[key] = []
			contentPropsAll[key].push({})
			// if (!contentPropsAll[key][j]) contentPropsAll[key][j] = {}
			if (val_str.match(/^http(s)?:\/\//)) {
				contentPropsAll[key][i].isUrl = true
				// row[key] = val.replace(/^http(s)?:\/\/([a-zA-Z0-9$-_.+!*'(),/&?=:%]+?)(\/)?$/, '$2') // Reformat URLs
			} else if (val_str.length > 70) {
				contentPropsAll[key][i].isLongText = true
			} else if (val === null) {
				contentPropsAll[key][i].empty = true
			}
		})
	})

	// Then for the sorters, we check if they're all the same,
	// and if they are, we set it as the column's sorter.
	// - - -
	// Supported sorter can be found here:
	// https://tabulator.info/docs/5.5/sort#func-builtin
	const sorters = {}
	Object.entries(sortersAll).forEach(([key, typesPerRow]) => {
		if (typesPerRow.every((type, i, arr) => type == arr[0])) {
			sorters[key] = typesPerRow[0]
		} else {
			sorters[key] = 'string'
		}
	})

	// Then for the content props, we check if at least
	// one cell in a column has props set, and if there
	// is, we set a custom formatter for that column,
	// which will format the data at the cell level,
	// based on the content of the individual cell.
	// - - -
	// The built-in formatters (eg. 'url') only work at
	// the column level, but we want more flexibility.
	// Eg. if one cell has long text, we need a textarea
	// editor for that cell, but we don't want to apply
	// this to all the column's cells.
	const needFormatter = {}
	Object.entries(contentPropsAll).forEach(([key, propsPerRow]) => {
		if (propsPerRow.some(props => !!Object.keys(props).length && !props.empty)) {
			needFormatter[key] = true
		} else {
			needFormatter[key] = false
		}
	})

	// There's a built-in 'textarea' editor, but we weren't happy with
	// how truncation was handled for cells containing long text, so
	// we had to create a custom formatter that gives us more control
	// (specifically that lets us truncate text to a custom number of
	// lines) and a custom formatter begets a custom editor.
	// - - -
	// https://tabulator.info/docs/5.5/edit#edit-custom
	function customTextareaEditor(cell, onRendered, success, cancel, editorParams) {
		const editor = document.createElement('textarea')
		editor.value = cell.getValue()

		onRendered(() => {
			// Expand textarea to fit content
			editor.style.height = ''
			editor.style.height = editor.scrollHeight + 'px'

			// Focus or select text
			if (editorParams.selectContents) {
				editor.select()
			} else {
				editor.focus()
			}
		})
		editor.addEventListener('change', () => {
			applyEdit()
		})
		editor.addEventListener('blur', () => {
			applyEdit()
		})
		document.addEventListener('keyup', e => {
			if (e.key == 'Enter') {
				applyEdit()
			} else if (e.key == 'Escape') {
				cancel()
			}
		})

		return editor

		//
		//

		function applyEdit() {
			success(editor.value)
		}
	}

	// Create columns
	const columns = Object.entries(sorters).map(([key, sorter], i) => {
		const { formatter, formatterParams, formatterType } = _pickFormatter(key, needFormatter, contentPropsAll)
		// About built-in editors: https://tabulator.info/docs/5.5/edit#edit-builtin
		const width = formatterType == 'textarea' ? 400 : undefined // Set column width only for columns with long text
		const editor = formatterType == 'textarea' ? customTextareaEditor : true // True will auto-select the editor, usually just 'input'
		const editorParams = { selectContents: true } // Select value on focus

		col = {
			sorter,
			title: key,
			field: key,
			width,
			editor,
			editorParams,
			editable: () => table.isEditMode(),
			formatter,
			formatterParams,
			headerMenu: contextMenus.header,
			headerContextMenu: contextMenus.header,

			// // This blocks the user from resizing the
			// // column, so we use a formatter instead:
			// maxWidth: 500,
		}

		// Add sorter params
		if (sorter == 'date') {
			col.sorterParams = {
				format: 'yyyy-MM-dd',
				alignEmptyValues: 'top',
			}
		} else if (sorter == 'number') {
			col.sorterParams = {
				thousandSeparator: ',',
				decimalSeparator: '.',
				alignEmptyValues: 'top',
			}
		}

		return col
	})

	console.log('Sorters:', sorters)
	console.log('Data:', data)
	console.log('Columns:', columns)
	console.log('ContentProps', contentPropsAll)
	console.log('needFormatter', needFormatter)
	return columns
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

// #endregion

/////////////////////////////////
// #region - Table interaction

// %% Trash
// table.on('headerClick', onHeaderClick)
// table.on('dataSorting', onDataSorting)

// function onHeaderClick(e, column) {
// 	e.preventDefault()
// 	console.log(111)
// }

// function onDataSorting(sorters) {
// 	console.log(222, sorters)
// 	return false
// }

function onCellClick(e, cell) {
	// Display overlay with full content when cell is truncated.
	if (e.target.classList.contains('text-wrap')) {
		const $textWrap = e.target
		const $cell = $textWrap.closest('.tabulator-cell')
		const isTruncated = _isTruncated($textWrap)
		// if (isTruncated) _expandCell($textWrap, $cell, cell)
		if (isTruncated) displayFullCellContent($textWrap, $cell)
	}

	// Set focus on clicked cell.
	$focusCell = cell.getElement()
	setFocus($focusCell)

	//
	//

	// Expand the cell so all content is visible (unused)
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
	if (e.shiftKey) {
		// Select all rows between the last selected row and the current row
		const selectedRows = table.getSelectedRows()
		if (selectedRows.length) {
			const lastSelectedRow = selectedRows[selectedRows.length - 1]
			const lastSelectedRowIndex = lastSelectedRow.getPosition()
			const currentRowIndex = row.getPosition()
			let lowIndex = Math.min(lastSelectedRowIndex, currentRowIndex)
			let highIndex = Math.max(lastSelectedRowIndex, currentRowIndex)

			// When you select from bottom to top, we gotta include the highIndex
			// When you select from top to bottom, we gotta include the lowIndex
			if (lowIndex != lastSelectedRowIndex) {
				lowIndex -= 1
				highIndex -= 1
			}

			const toBeSelected = table.getRows().slice(lowIndex, highIndex)
			console.log(toBeSelected.map(row => row.getPosition()))
			if (table.lastClickedRowSelected) {
				console.log('A')
				table.selectRow(toBeSelected)
				table.lastClickedRowSelected = true
			} else {
				console.log('B')
				table.deselectRow(toBeSelected)
				table.lastClickedRowSelected = false
			}
		}
	} else {
		row.toggleSelect() //toggle row selected state on row click
		table.lastClickedRowSelected = table.getSelectedRows().includes(row)
	}
	console.log('###', table.lastClickedRowSelected)
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
		if (e.relatedTarget && e.relatedTarget.getElementById('display-full-text')) {
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
function deselectRows() {
	const selectedRows = table.getSelectedRows()
	if (selectedRows.length > 3) {
		if (confirm('Are you sure you want to delect all rows?')) {
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

// Copies the content of a cell to the clipboard
function copyCellContent($cell) {
	if (navigator.clipboard) {
		// New way
		navigator.clipboard.writeText($cell.innerText)
	} else {
		// Old way
		const range = document.createRange()
		range.selectNode($cell)
		window.getSelection().removeAllRanges()
		window.getSelection().addRange(range)
		document.execCommand('copy')
		window.getSelection().removeAllRanges()
	}

	$cell.classList.add('copied')
	setTimeout(() => {
		$cell.classList.remove('copied')
	}, 600)
}

// %% trash
// // Generate the interactive HTML for the context menu
// // Note: there's a built-in contextMenu / headerContextMenu
// // which we could have used instead, but is has some issues
// function showContextMenu(col) {
// 	const colName = col.getField()
// 	$headerCol = document.querySelector(`#table.tabulator .tabulator-header .tabulator-col[tabulator-field="${colName}"]`)

// 	// Remove any other context menu
// 	const $prevMenu = document.getElementById('context-menu')
// 	if ($prevMenu) {
// 		$prevMenu.remove()
// 	}

// 	// Create HTML
// 	const $contextMenu = document.createElement('div')
// 	$contextMenu.setAttribute('id', 'context-menu')
// 	$contextMenu.innerHTML = `<div class="menu-item" value="col-remove" data-col-name="${colName}">Remove column</div>`

// 	// Position right below the clicked header column
// 	const rect = $headerCol.getBoundingClientRect()
// 	$contextMenu.style.setProperty('top', rect.bottom + 'px')
// 	$contextMenu.style.setProperty('left', rect.left + 'px')
// 	$contextMenu.style.setProperty('width', rect.width - 1 + 'px')

// 	document.body.append($contextMenu)
// }

// // Context menu click handler
// function onContextMenuClick(e) {
// 	if (e.target.classList.contains('menu-item')) {
// 		if (e.target.getAttribute('value') == 'col-remove') {
// 			const colName = e.target.getAttribute('data-col-name')
// 			table.deleteColumn(colName)
// 		}
// 	}

// 	// Remove context menu from DOM
// 	const $contextMenu = document.getElementById('context-menu')
// 	if ($contextMenu) {
// 		$contextMenu.remove()
// 	}
// }

// %% trash
// This is old code for column context menu, but may be useful
// if we want to add other context menus.
// function initRightClickMenu() {
// 	// Add event listener to each column header
// 	const $headerCols = document.querySelectorAll('#table.tabulator .tabulator-header .tabulator-col')
// 	$headerCols.forEach($headerCol => {
// 		$headerCol.addEventListener('contextmenu', function (e) {
// 			_renderContextMenu($headerCol)
// 			e.preventDefault()
// 		})
// 	})

// 	//
// 	//

// 	// Generate the interactive HTML for the context menu
// 	function _renderContextMenu($headerCol) {
// 		const colName = $headerCol.querySelector('.tabulator-col-title').innerText

// 		// Remove any other context menu
// 		const $prevMenu = document.getElementById('context-menu')
// 		if ($prevMenu) {
// 			$prevMenu.remove()
// 		}

// 		// Create HTML
// 		const $menu = document.createElement('div')
// 		$menu.setAttribute('id', 'context-menu')
// 		$menu.innerHTML = '<div class="menu-item" value="remove">Remove column</div>'

// 		// Click handler
// 		$menu.addEventListener('click', e => {
// 			if (e.target.classList.contains('menu-item')) {
// 				if (e.target.getAttribute('value') == 'remove') {
// 					console.log('Remove column', colName)
// 					table.deleteColumn(colName)
// 				} else if (e.target.getAttribute('value') == 'reset') {
// 					console.log('Reset columns')
// 				} else if (e.target.getAttribute('value') == 'foo') {
// 					console.log('Foo columns')
// 				}
// 			}
// 			$menu.remove()
// 		})

// 		// Position right below the clicked header column
// 		const rect = $headerCol.getBoundingClientRect()
// 		$menu.style.setProperty('top', rect.bottom + 'px')
// 		$menu.style.setProperty('left', rect.left + 'px')
// 		$menu.style.setProperty('width', rect.width - 1 + 'px')

// 		document.body.append($menu)

// 		return $menu
// 	}
// }

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

// Set truncation limit (0/1/3)
function setTruncation(limit) {
	const $table = document.getElementById('table')
	$table.classList.remove('trunc-3', 'trunc-1')
	if (limit > 0) {
		$table.classList.add('trunc-' + limit)
	}
	table.redraw(true)
}

// Apply options
function applyOptions() {
	// Truncation
	const truncLimit = +document.getElementById('opt-select-truncation').value
	setTruncation(truncLimit)

	// Index column
	const includeIndexCol = +document.getElementById('opt-add-index-col').value
	if (includeIndexCol) {
		table.addIndexCol(true)
	} else {
		table.removeIndexCol()
	}
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

// Set the "Add index column" dropdown to the correct value.
function setIndexColDropdown() {
	if (table.addedIndex) {
		document.getElementById('opt-add-index-col').querySelector('option[value="1"]').selected = true
	} else {
		document.getElementById('opt-add-index-col').querySelector('option[value="0"]').selected = true
	}
}

// #endregion
