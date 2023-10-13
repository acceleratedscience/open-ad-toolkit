/**
 * Click handlers
 */

// Edit
document.getElementById('btn-edit').addEventListener('click', function () {
	table.toggleEditMode(true)
})

// Cancel
document.querySelector('#btn-cancel').addEventListener('click', function () {
	table.toggleEditMode(false)
})

// Save
document.querySelector('#btn-save').addEventListener('click', function () {
	table.toggleEditMode(false)
})

// Reset columns
document.querySelector('#table-wrap > a').addEventListener('click', function () {
	table.resetCols()
})

// Move focus with arrow keys
document.addEventListener('keydown', e => {
	if (!document.querySelector('.tabulator-cell.focus')) return
	if (isEditMode()) return

	if (e.key == 'ArrowLeft') {
		moveFocus('left')
	} else if (e.key == 'ArrowRight') {
		moveFocus('right')
	} else if (e.key == 'ArrowUp') {
		moveFocus('up')
	} else if (e.key == 'ArrowDown') {
		moveFocus('down')
	} else if (e.key == 'Escape') {
		document.querySelector('.tabulator-cell.focus').classList.remove('focus')
	} else if (e.key == 'Enter') {
		table.toggleEditMode(true)
		setTimeout(() => {
			document.querySelector('.tabulator-cell.focus').click()
		}, 100)
	}
	e.preventDefault()
})

// Exit focus on blur
document.addEventListener('click', e => {
	if (e.target.closest('.tabulator-cell') || e.target.classList.contains('tabulator-cell')) return
	const $currentFocusCell = document.querySelector('.tabulator-cell.focus')
	if ($currentFocusCell) $currentFocusCell.classList.remove('focus')
})

// Copy focused cell content on cmd + C
document.addEventListener('keydown', e => {
	if (e.key == 'c' && (e.metaKey || e.ctrlKey)) {
		$cell = document.querySelector('.tabulator-cell.focus')
		copyCellContent($cell)
	}
})

// Truncation setting
document.getElementById('select-truncation').addEventListener('change', setTruncation)

/**
 * Table rendering
 */

// Parse data
const data = JSON.parse(document.getElementById('table').getAttribute('data'))
document.getElementById('table').removeAttribute('data')

// Add index column if not yet included
data.forEach((row, i) => {
	let hasIndex = !!row['#']
	if (!hasIndex) {
		for (const key in row) {
			if (key.toLowerCase() == 'index') {
				hasIndex = true
				break
			}
		}
	}
	if (!hasIndex) {
		row['#'] = i + 1
	}
})

// Parse columns
const columns = parseColumns(data)
// columns.unshift({ formatter: 'rowSelection', titleFormatter: 'rowSelection', align: 'center', headerSort: false }) // Add checkbox column fopr selection UI

// Create table
const table = new Table('#table', {
	data,
	columns,
	resizableRows: true,
	// rowHeight: 40,
	reactiveData: true,
	editMode: window.location.hash == '#edit',
	pagination: false,
	movableRows: true,
	movableColumns: true,
	// selectable: true,
	// paginationSize: 5,
	// rowFormatter: row => {
	// 	// Trash
	// 	// row.getColumns()
	// 	// const rowData = row.getData()
	// 	// Calculate the row height based on your criteria
	// 	// const height = 26
	// 	// Set the row height using CSS
	// 	// row.getElement().style.height = height + 'px'
	// },
})

// Redraw the table after the data object is built but before the table is rendered.
// This is a little hack to prevent all rows to take on the height of the tallest cell.
table.on('tableBuilt', () => {
	table.redraw(true)
})
table.on('cellClick', onCellClick)
table.on('cellDblClick', e => {
	$cell = e.target.classList.contains('tabulator-cell') ? e.target : e.target.closest('.tabulator-cell')
	copyCellContent($cell)
})

//
//
//
//

// Parse data & create columns.
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

	// %%
	// https://tabulator.info/docs/5.5/edit#edit-custom
	var customEditor = function (cell, onRendered, success, cancel, editorParams) {
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
		const width = formatterType == 'textarea' ? 400 : undefined
		const editor = formatterType == 'textarea' ? customEditor : true
		const editorParams = { selectContents: true }

		col = {
			sorter,
			title: key,
			field: key,
			editor,
			editorParams,
			editable: isEditMode,
			// // This blocks the user from resizing the
			// // column, so we use a formatter instead:
			// maxWidth: 500,

			formatter,
			formatterParams,
			width,
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

	// Find the column with field name 'index' (cas insensitive) or '#'
	// and move it to the front of the array.
	const indexColIndex = columns.findIndex(col => {
		return col.field.toLowerCase() == 'index' || col.field == '#'
	})
	const indexCol = columns.splice(indexColIndex, 1)[0]
	columns.unshift(indexCol)

	console.log('Sorters:', sorters)
	console.log('Data:', data)
	console.log('Columns:', columns)
	console.log('ContentProps', contentPropsAll)
	console.log('needFormatter', needFormatter)
	return columns
}

function _pickFormatter(key, needFormatter, contentPropsAll) {
	hasSomeLongText = contentPropsAll[key].some(props => props.isLongText)
	isAllUrls = contentPropsAll[key].every(props => props.isUrl)

	if (needFormatter[key]) {
		if (isAllUrls) {
			return { formatter: 'link', formatterParams: { target: '_blank' } }
		} else if (hasSomeLongText) {
			// return { formatter: 'textarea' } // %%
			return {
				formatterType: 'textarea',
				formatter: (cell, formatterParams, onRendered) => {
					return `<div class="text-wrap">${cell.getValue()}</div>`
				},
			}
		}
	}
	return {}
}

function onCellClick(e, cell) {
	// Display overlay with full content when cell is truncated.
	if (e.target.classList.contains('text-wrap')) {
		const $textWrap = e.target
		const $cell = $textWrap.closest('.tabulator-cell')
		const isTruncated = _isTruncated($textWrap)
		// if (isTruncated) _expandCell($textWrap, $cell, cell)
		if (isTruncated) _displayFullContent($textWrap, $cell)
	}

	// Set focus on clicked cell.
	$focusCell = cell._cell.element
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

	// Display overlay div that matches the cell's content and position,
	// but that expands to the bottom to display the full, untruncated content.
	function _displayFullContent($textWrap, $cell, $row, row) {
		const display = document.createElement('div')
		display.classList.add('display-full-text')
		display.setAttribute('tabindex', 0)
		display.innerHTML = $textWrap.innerHTML
		$cell.append(display)
		display.focus()
		display.addEventListener('blur', _removeDisplay)
	}

	function _removeDisplay(e) {
		// Don't remove display when display itself is clicked.
		if (e.relatedTarget && e.relatedTarget.querySelector('.display-full-text')) {
			e.target.focus()
			return
		}

		e.target.remove()
		delete e.target
	}
}

// Set artificial cell focus.
function setFocus($focusCell) {
	const $currentFocusCell = document.querySelector('.tabulator-cell.focus')
	if ($currentFocusCell) $currentFocusCell.classList.remove('focus')
	$focusCell.classList.add('focus')
}

// Move artificial cell focus.
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

/**
 * Utility functions
 */

// Check if string is a date.
function isDate(str, log) {
	// Ignore numbers that are not formatted like a date.
	const separators = str.match(/[-./]/g, '')
	if (log) console.log(separators)
	if (separators && separators.length == 2 && separators[0] == separators[1]) {
		return moment(str).isValid()
	} else if (str.match(/[jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec]/i)) {
		return moment(str).isValid()
	} else {
		return false
	}
}

// Check if table is in edit mode.
function isEditMode(cell, a) {
	return table.editMode
}

// Set truncation limit (0/1/3)
function setTruncation(e) {
	truncLimit = e.target.value
	const $table = document.getElementById('table')
	$table.classList.remove('trunc-3', 'trunc-1')
	if (+truncLimit > 0) {
		$table.classList.add('trunc-' + truncLimit)
	}
	table.redraw(true)
}
