/**
 * Table is a subclass that extends
 * Tabulator with some extra functionality.
 */

class Table extends Tabulator {
	constructor(element, options) {
		// Initiation
		super(element, options)
		super.on('tableBuilt', () => {
			// Turn on edit mode if hash is present
			if (options.editMode) {
				this.toggleEditMode(true)
			}

			// Store column widths so we can reset them
			this.storeColWidths()
		})

		// Modifation
		super.on('columnResized', () => {
			this.element.classList.add('resized-col')
			table.redraw(true)
		})
		// super.on('rowResized', () => { // %% trash
		// 	this.element.classList.add('resized-row')
		// })
		super.on('dataProcessed', () => {
			if (!this.initialized) {
				this.originalData = this.getData() // See hasEdits()
				this.initialized = true
			}
		})

		// Variables
		this.initialized = false // Has the table been initialized
		this.originalData = null // See hasEdits()
		this.editMode = false // Used to determine if we're in edit mode - see isEditMode()
		this.colDefaultWidths = {} // Used to store column widths so we can reset them
		// this.rowDefaultHeights = {} // Used to store row heights so we can reset them
		this.index = options.index // Tabulator doesn't store the index so we have to
		this.addedIndex = options.addedIndex // Used to keep track if we added an index column
		this.lastSelectedRowIndex = null // Number used to determine where to start from when shift-selecting
		this.lastSelectedRowSelState = null // Boolean used to determine if shift-select row should batch-select or batch-deselect
		this.selectMode = false // Boolean used to ignore row resize handles while selecting rows
		this.history = [] // Used to store previous 10 states so you can undo/redo
		this.history_reverse = [] // Used to store the states you went back in history
	}

	// Store column default widths
	storeColWidths() {
		this.getColumns().forEach(col => {
			this.colDefaultWidths[col.getField()] = col.getWidth()
		})
	}

	// // Store row default heights
	// storeRowHeights() {
	// 	this.getRows().forEach(row => {
	// 		this.rowDefaultHeights[row.getPosition()] = row.getHeight()
	// 	})
	// }

	// // Add index column
	// addIndexCol(force) {
	// 	if (this.addedIndex) return
	// 	const data = this.getData()

	// 	let addCol = false
	// 	if (force) {
	// 		// Column is added manually via options panel.
	// 		addCol = true
	// 	} else if (!force) {
	// 		// Column is added on init, but only if one doesn't already exist.
	// 		const rowSample = data[1]
	// 		let hasIndex = !!rowSample['#']
	// 		if (!hasIndex) {
	// 			for (const key in rowSample) {
	// 				if (key.toLowerCase() == 'index') {
	// 					hasIndex = true
	// 					break
	// 				}
	// 			}
	// 			addCol = !hasIndex
	// 		}
	// 	}

	// 	// Add index column if missing
	// 	if (addCol) {
	// 		this.addedIndex = true
	// 		data.forEach((row, i) => {
	// 			row['#'] = i + 1
	// 		})
	// 		this.setData(data)

	// 		const colSample = this.getColumns()[0]
	// 		const newCol = {
	// 			...colSample,
	// 			sorter: 'number',
	// 			title: '#',
	// 			field: '#',
	// 		}
	// 		this.addColumn(newCol, true)
	// 	}
	// }

	// // Remove index column
	// removeIndexCol() {
	// 	const indexCol = this.getColumns()[0]
	// 	if (indexCol.getField() == '#') {
	// 		this.addedIndex = false
	// 		indexCol.delete()
	// 		const data = this.getData()
	// 		data.forEach((row, i) => {
	// 			delete row['#']
	// 		})
	// 		this.setData(data)
	// 	}
	// }

	// Reset column widths to default
	resetCols() {
		this.getColumns().forEach(col => {
			col.setWidth(this.colDefaultWidths[col.getField()])
			this.element.classList.remove('resized-col')
		})
	}

	// // Reset row heights to default
	// resetRows() {
	// 	this.getRows().forEach(row => {
	// 		row.setHeight(this.rowDefaultHeights[row.getField()])
	// 		this.element.classList.remove('resized-row')
	// 	})
	// }

	// Check if table is in edit mode
	isEditMode() {
		return this.editMode
	}

	// Toggle edit mode
	// - - -
	// Note: Tabulator doesn't support toggling edit mode, so we use a little hack.
	// The actual toggling of edit mode is done via isEditMode when creating the table.
	toggleEditMode(bool, revertChanges) {
		this.editMode = bool == undefined ? !this.editMode : bool
		if (this.editMode) {
			// ENTER
			this.storeData() // Store data so we can revert on cancel
			this.element.classList.add('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search + '#edit') // Add hash
		} else {
			// EXIT
			if (revertChanges) this.revertData()
			this.element.classList.remove('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search) // Remove hash
		}
	}

	// Store copy of data so we can revert changes
	storeData() {
		this.dataBeforeEdit = this.getData()
	}

	// Revert changes on cancel
	async revertData() {
		if (this.hasEdits(true)) {
			const selectedRows = await this.getSelectedRows()
			const selectedRowsIndexes = selectedRows.map(row => row.getIndex())
			this.setData(this.dataBeforeEdit)
			this.redraw()

			// Re-select rows that were selected before
			const newRows = table.getRows()
			const newSelectedRows = selectedRowsIndexes.map(index => newRows[index - 1])
			this.selectRow(newSelectedRows)
		}
	}

	// Check if table has edits:
	// - In this session only (i.e. since last save) --> used to revert data on cancel
	// - Since the table was loaded --> used to prevent accidental page exit
	hasEdits(sessionOnly) {
		const dataBefore = sessionOnly ? this.dataBeforeEdit : this.originalData
		const dataEdited = this.getData()
		if (dataBefore.length != dataEdited.length) return true
		for (let i = 0; i < dataBefore.length; i++) {
			const row1 = dataBefore[i]
			const row2 = dataEdited[i]
			for (const key in row1) {
				if (row1[key] != row2[key]) return true
			}
		}
		return false
	}

	// // %% trash
	// // Recreate the built in sort
	// sort(col) {
	// 	// We block sorting while in edit mode, because we reset the selection
	// 	// on revertData() and sorting breaks this. However, this is not a good
	// 	// solution, instead we should detect if the data was sorted, and not
	// 	// reselect the rows if it was.
	// 	// if (this.editMode) {
	// 	// 	alert('Please exit edit mode before sorting.')
	// 	// 	return
	// 	// }
	// 	const colName = col.getField()
	// 	const $col = col.getElement()
	// 	const sortOrder = $col.getAttribute('aria-sort') == 'ascending' ? 'desc' : 'asc'
	// 	$col.classList.add('sortable', `sort-${sortOrder}`)
	// 	this.setSort(colName, sortOrder)
	// 	console.log('$')
	// }

	// Homebrewed version of movableRows.
	// - - -
	// We had to build our own because the build-in movableRows
	// can't be made conditional, and it creates undesireable
	// side effects in edit mode. Specifically, it will interfere
	// with the textarea resize handle.
	onCellMouseDown(clickEvent, cell, onCellClick) {
		if (this.editMode) return

		// Only accept left click.
		if (clickEvent.button > 0) return

		// This lets us check if it's a drag or a click
		let dragging = false
		const debug = false // Display indicators for debugging

		// Store references
		const table = this
		const row = cell.getRow()
		const $row = row.getElement()
		const $rowClone = $row.cloneNode(true)
		const $parent = this.element.querySelector('.tabulator-tableholder')

		// Calculate values
		const rowY = $row.getBoundingClientRect().top
		let parentY = $parent.getBoundingClientRect().top
		let y = rowY - parentY - 1
		let mouseRowOffset = clickEvent.clientY - rowY // Offset between cursor and row top
		let jump = 0 // How many positions the row has jumped

		// Listen to mouse movement
		// - - -
		// If mouse position moves before the mouseup event fires,
		// we start moving the row. Otherwise it registers as a click.
		document.addEventListener('mousemove', _onMouseMove)
		document.addEventListener('mouseup', _onMouseUp)

		// Debug
		let $reference1
		let $reference2

		// Drag listener
		function _onMouseMove(e) {
			if (!dragging) {
				dragging = true
				_onMouseMoveInit()
			}
			parentY = $parent.getBoundingClientRect().top
			y = e.clientY - parentY - mouseRowOffset
			const minY = -1
			const maxY = $parent.clientHeight - $rowClone.clientHeight - 2
			y = Math.max(Math.min(y, maxY), minY)
			$rowClone.style.top = y + 'px'
			_moveRow(y)
		}

		// Executes once when you start dragging
		function _onMouseMoveInit() {
			// Update DOM
			$row.classList.add('while-dragging')
			$rowClone.classList.remove('tabulator-selectable')
			$rowClone.classList.add('dragging')
			$rowClone.style.top = y + 'px'
			$parent.appendChild($rowClone)

			// Debug
			if (debug) {
				$reference1 = document.createElement('div')
				$reference1.style.height = '1px'
				$reference1.style.width = '100px'
				$reference1.style.background = 'blue'
				$reference1.style.position = 'absolute'
				$reference1.style.left = 0
				$reference1.style.zIndex = 10
				$parent.appendChild($reference1)
				//
				$reference2 = document.createElement('div')
				$reference2.style.height = '1px'
				$reference2.style.width = '100px'
				$reference2.style.background = 'green'
				$reference2.style.position = 'absolute'
				$reference2.style.left = 0
				$reference2.style.zIndex = 10
				$parent.appendChild($reference2)
			}
		}

		// On drag - move row to closest position
		function _moveRow(y) {
			const yRowMiddle = y + $rowClone.clientHeight / 2
			const moveRow = _findTargetIndex(yRowMiddle)
			jump += moveRow

			if (moveRow == 1) {
				const $nextRow = $row.nextElementSibling
				$nextRow.after($row)
			} else if (moveRow == -1) {
				const $prevRow = $row.previousElementSibling
				$prevRow.before($row)
			}
		}

		// Return -1 / 0 / 1 depending on where the row should be moved
		function _findTargetIndex(yRowMiddle) {
			const $prevRow = $row.previousElementSibling
			const $nextRow = $row.nextElementSibling
			const prevRowMiddle = $prevRow ? $prevRow.offsetTop + $prevRow.clientHeight / 2 : null
			const nextRowMiddle = $nextRow ? $nextRow.offsetTop + $nextRow.clientHeight / 2 : null
			const thisRowTop = $rowClone.offsetTop
			const thisRowBottom = $rowClone.offsetTop + $row.clientHeight

			// Debug
			if (debug) {
				$reference1.style.top = prevRowMiddle + 'px'
				$reference2.style.top = nextRowMiddle + 'px'
			}

			if (prevRowMiddle != null && thisRowTop < prevRowMiddle) {
				return -1
			} else if (nextRowMiddle != null && thisRowBottom > nextRowMiddle) {
				return 1
			} else {
				return 0
			}
		}

		// On drop
		function _onMouseUp() {
			document.removeEventListener('mousemove', _onMouseMove)
			document.removeEventListener('mouseup', _onMouseUp)

			// If no dragging occured, it's a click.
			if (!dragging) {
				onCellClick(clickEvent, cell)
				return
			}

			$rowClone.remove()
			$row.classList.remove('while-dragging')

			const data = table.getData()
			const currentPos = row.getPosition() - 1
			const newPos = currentPos + jump
			data.splice(newPos, 0, data.splice(currentPos, 1)[0])
			table.updateData(data)

			// Debug
			if (debug && $reference1 && $reference2) {
				console.log(currentPos, '-->', newPos)
				console.log(data.map(itm => itm.Index))
				console.log(data)
				$reference1.remove()
				$reference2.remove()
			}
		}
	}

	// moveRowStart(e, row) {
	// 	if (this.editMode) return

	// 	// Store elements
	// 	const $row = row.getElement()
	// 	const $rowClone = $row.cloneNode(true)
	// 	const $parent = this.element.querySelector('.tabulator-tableholder')

	// 	// Calculate values
	// 	const rowY = $row.getBoundingClientRect().top
	// 	let parentY = $parent.getBoundingClientRect().top
	// 	let y = rowY - parentY - 1
	// 	const mouseRowOffset = e.clientY - rowY
	// 	let jump = 0 // How many positions the row has jumped
	// 	const table = this

	// 	// Update DOM
	// 	$row.classList.add('while-dragging')
	// 	$rowClone.classList.remove('tabulator-selectable')
	// 	$rowClone.classList.add('dragging')
	// 	$rowClone.style.top = y + 'px'
	// 	$parent.appendChild($rowClone)

	// 	document.addEventListener('mousemove', _moveRowDrag)
	// 	document.addEventListener('mouseup', _moveRowEnd)

	// 	// // For testing
	// 	// const $indicator = document.createElement('div')
	// 	// $indicator.style.height = '1px'
	// 	// $indicator.style.width = '100px'
	// 	// $indicator.style.background = 'red'
	// 	// $indicator.style.position = 'absolute'
	// 	// $indicator.style.zIndex = 1000
	// 	// $parent.appendChild($indicator)
	// 	// const $reference1 = document.createElement('div')
	// 	// $reference1.style.height = '1px'
	// 	// $reference1.style.width = '100px'
	// 	// $reference1.style.background = 'blue'
	// 	// $reference1.style.position = 'absolute'
	// 	// $reference1.style.left = 0
	// 	// $reference1.style.zIndex = 1000
	// 	// $parent.appendChild($reference1)
	// 	// const $reference2 = document.createElement('div')
	// 	// $reference2.style.height = '1px'
	// 	// $reference2.style.width = '100px'
	// 	// $reference2.style.background = 'green'
	// 	// $reference2.style.position = 'absolute'
	// 	// $reference2.style.left = 0
	// 	// $reference2.style.zIndex = 1000
	// 	// $parent.appendChild($reference2)

	// 	// On drag listener
	// 	function _moveRowDrag(e) {
	// 		parentY = $parent.getBoundingClientRect().top
	// 		y = e.clientY - parentY - mouseRowOffset
	// 		const minY = -1
	// 		const maxY = $parent.clientHeight - $rowClone.clientHeight - 2
	// 		y = Math.max(Math.min(y, maxY), minY)
	// 		$rowClone.style.top = y + 'px'
	// 		_moveRow(y)
	// 	}

	// 	// On drag - move row to closest position
	// 	function _moveRow(y) {
	// 		const yRowMiddle = y + $rowClone.clientHeight / 2
	// 		const moveRow = _findTargetIndex(yRowMiddle)
	// 		jump += moveRow

	// 		if (moveRow == 1) {
	// 			const $nextRow = $row.nextElementSibling
	// 			$nextRow.after($row)
	// 		} else if (moveRow == -1) {
	// 			const $prevRow = $row.previousElementSibling
	// 			$prevRow.before($row)
	// 		}
	// 	}

	// 	function _findTargetIndex(yRowMiddle) {
	// 		const $prevRow = $row.previousElementSibling
	// 		const $nextRow = $row.nextElementSibling
	// 		const prevRowMiddle = $prevRow ? $prevRow.offsetTop + $prevRow.clientHeight / 2 : null
	// 		const nextRowMiddle = $nextRow ? $nextRow.offsetTop + $nextRow.clientHeight / 2 : null
	// 		const thisRowTop = $rowClone.offsetTop
	// 		const thisRowBottom = $rowClone.offsetTop + $row.clientHeight

	// 		// // For testing
	// 		// $indicator.style.top = yRowMiddle + 'px'
	// 		// $reference1.style.top = prevRowMiddle + 'px'
	// 		// $reference2.style.top = nextRowMiddle + 'px'

	// 		if (prevRowMiddle != null && thisRowTop < prevRowMiddle) {
	// 			return -1
	// 		} else if (nextRowMiddle != null && thisRowBottom > nextRowMiddle) {
	// 			return 1
	// 		} else {
	// 			return 0
	// 		}
	// 	}

	// 	// On drop
	// 	function _moveRowEnd() {
	// 		document.removeEventListener('mousemove', _moveRowDrag)
	// 		document.removeEventListener('mouseup', _moveRowEnd)
	// 		$rowClone.remove()
	// 		$row.classList.remove('while-dragging')

	// 		const data = table.getData()
	// 		// console.log('')
	// 		// console.log(data.map(itm => itm.Index))
	// 		const currentPos = row.getPosition() - 1
	// 		const newPos = currentPos + jump
	// 		data.splice(newPos, 0, data.splice(currentPos, 1)[0])
	// 		table.setData(data)
	// 		// console.log(currentPos, '-->', newPos)
	// 		// console.log(data.map(itm => itm.Index))

	// 		// // For testing
	// 		// $indicator.remove()
	// 		// $reference1.remove()
	// 		// $reference2.remove()
	// 	}
	// }

	// Returns an object with arrays of cell property objects per row:
	// { col1: [
	// 	  { type: 'string', isUrl: false, isLongText: false }
	// 	  { type: 'string', isUrl: true, isLongText: false }
	// 	  { type: 'string', isUrl: false, isLongText: true }
	//   ],
	//   col2: ...
	// }
	getCellProps(data) {
		data = data ? data : table.getData()
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

	// Add state of data to history
	add_history(data) {
		data = JSON.parse(JSON.stringify(data))
		console.log(111, data['Organization Id:'], data)
		if (this.history.length >= 10) {
			this.history.shift()
		}
		this.history.push(data)
	}

	undo() {
		if (this.history.length > 0) {
			const dataCurrent = this.history.pop()
			this.history_reverse.push(dataCurrent)
			const dataPrevious = this.history.length > 0 ? this.history[this.history.length - 1] : null
			this.setData(dataPrevious)
		}
	}

	redo() {
		if (this.history_reverse.length > 0) {
			const data = this.history_reverse.pop()
			this.history.push(data)
			this.updateData(data)
		}
	}

	// Tag on to the Tabulator redraw function.
	redraw(bool) {
		console.log('redraw', bool)
		if (!bool) {
			super.redraw()
			return
		}

		// Store current column widths.
		const colWidthsBefore = {}
		table.getColumns().forEach(col => {
			colWidthsBefore[col.getField()] = col.getWidth()
		})

		// Redraw
		super.redraw(bool)

		// Store new column widths.
		const colWidthsAfter = {}
		table.getColumns().forEach(col => {
			colWidthsAfter[col.getField()] = col.getWidth()
		})

		// Fix width per column.
		this._fixColumnWidths(colWidthsBefore, colWidthsAfter)
	}

	// Fix the width per column, so the table doesn't snap back
	// to the automatic content-based width and truncation is applied.
	// If the column has not been resized manually, we'll limit
	// it to 400px, otherwise we'll keep the width as it was.
	// - - -
	// Note: Tabulator has a built-in max-width property, but we
	// don't want to use this because it prevents the user from
	// resizing a column past this width.
	_fixColumnWidths(colWidthsBefore, colWidthsAfter) {
		Object.entries(colWidthsAfter).forEach(([colName, width], i) => {
			const isUntampered = colWidthsBefore[colName] == this.colDefaultWidths[colName]
			// prettier-ignore
			const newWidth = isUntampered
				? (Math.min(colWidthsAfter[colName], 400))
				: (colWidthsBefore[colName])
			this.getColumn(colName).setWidth(newWidth)
		})
	}
}

Table.moduleName = 'custom'
