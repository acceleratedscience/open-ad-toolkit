/**
 * Table is a subclass that extends
 * Tabulator with some extra functionality.
 */

class Table extends Tabulator {
	constructor(element, options) {
		// Initiate Tabulator
		super(element, options)

		/**
		 * Events
		 */

		// Initialization
		super.on('tableBuilt', () => {
			// Redraw the table after the data object is built but before the table is rendered.
			// This is a little hack to prevent all rows to take on the height of the tallest cell.
			this.redraw(true)

			// Turn on edit mode if hash is present
			if (options.editMode) {
				this.toggleEditMode(true)
			}

			// Store column widths so we can reset them
			this._storeColWidths()

			// Store the initial state of the table.
			this._addHistoryEntry(this.getData())

			// Used to check if the table has been edited - See hasEdits()
			this.originalData = this.getData()

			// Keyboard interaction
			const onKeyDown = this.onKeyDown.bind(this)
			document.removeEventListener('keydown', onKeyDown)
			document.addEventListener('keydown', onKeyDown)
		})

		// Undo/redo
		super.on('dataChanged', data => {
			if (!this.blockHistory) {
				this._addHistoryEntry(data)
			}
		})

		// Column resizing
		super.on('columnResized', col => {
			const field = col.getField()
			if (!this.tamperedCols.includes(field)) {
				this.tamperedCols.push(field)
			}
			// Show "Reset columns" link
			this.element.classList.add('resized-col')

			// Updates row height to fit content.
			this.redraw(true)
		})

		// Controls our homebrewed movableRows.
		super.on('cellMouseDown', this.onCellMouseDown)

		// Update truncation whenever text is changed.
		super.on('cellEdited', cell => {
			const index = cell.getRow().getIndex()
			const field = cell.getField()

			// This will adjust the column widths to fit the content.
			this.redraw(true)

			// Trace back the same cell after redraw
			// has replaced the HTML and reset the focus.
			const $newCell = this.getRows()[index - 1].getCell(field).getElement()
			this.setFocus($newCell)
		})

		// Double-click to edit cell.
		super.on('cellDblClick', (e, cell) => {
			if (!this.isEditMode()) toggleEditMode(true)
			const $focusCell = cell.getElement()
			setTimeout(() => {
				$focusCell.click()
			}, 0)
		})

		/**
		 * Variables
		 */

		this.originalData = null // See hasEdits()
		this.editMode = false // Used to determine if we're in edit mode - see isEditMode()
		this.colDefaultWidths = {} // Used to store column widths so we can reset them
		// this.index = options.index // Tabulator doesn't store the index so we have to // %% trash --> this.options.index
		this.addedIndex = options.addedIndex // Used to keep track if we added an index column
		this.lastSelectedRowIndex = null // Number used to determine where to start from when shift-selecting
		this.lastSelectedRowSelState = null // Boolean used to determine if shift-select row should batch-select or batch-deselect
		this.selectMode = false // Boolean used to ignore row resize handles while selecting rows
		this.history = [] // Used to store previous 10 states so you can undo/redo
		this.history_reverse = [] // Used to store the states you went back in history
		this.blockHistory = false // Used to block writing to history when undoing/redoing
		this.isSubmitting = false // Used to prevent triggering window.onbeforeunload
		this.tamperedCols = [] // To keep track which columns have been manually resized
	}

	/////////////////////////////////
	// #region - Initialization

	// Store column default widths
	_storeColWidths() {
		this.getColumns().forEach(col => {
			this.colDefaultWidths[col.getField()] = col.getWidth()
		})
	}

	// #endregion

	/////////////////////////////////
	// #region - Public methods

	/**
	 * Editing
	 */

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
			this._storeData() // Store data so we can revert on cancel
			this.element.classList.add('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search + '#edit') // Add hash
		} else {
			// EXIT
			if (revertChanges) this._revertData()
			this.element.classList.remove('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search) // Remove hash
		}
	}

	// Store copy of data so we can revert changes
	_storeData() {
		this.dataBeforeEdit = this.getData()
	}

	// Revert changes on cancel
	async _revertData() {
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

	/**
	 * Column actions
	 */

	// Reset column widths to default
	resetColWidths() {
		this.getColumns().forEach(col => {
			col.setWidth(this.colDefaultWidths[col.getField()])
			this.element.classList.remove('resized-col')
		})
	}

	// Delete column and all its data
	deleteColumn(e, col) {
		col.delete()
	}

	// Rename column
	renameColumn(e, col) {
		const currentName = col.getField()
		const newName = prompt('New column name:', currentName)
		const cellWidth = col.getWidth()
		if (newName) {
			table.updateColumnDefinition(currentName, { title: newName, field: newName })
			table.getColumn(newName).setWidth(cellWidth)

			// Update every row with the new column name.
			table.getRows().forEach(row => {
				const updateObject = {}
				updateObject[currentName] = undefined
				updateObject[newName] = row.getCell(currentName).getValue()
				row.update(updateObject)
			})
		}
	}

	/**
	 * Undo / Redo
	 */

	// Add state of data to history
	_addHistoryEntry(data) {
		console.log('_addHistoryEntry')
		data = JSON.parse(JSON.stringify(data))
		if (this.history.length >= 10) {
			this.history.shift()
		}
		this.history.push(data)

		// Empty history_reverse so you can't redo
		// from another branch of edits.
		this.history_reverse = []
	}

	// Undo action
	async undoChange() {
		if (this.history.length > 1) {
			const dataCurrent = this.history.pop()
			this.history_reverse.push(dataCurrent)
			const dataPrevious = this.history[this.history.length - 1]
			this.blockHistory = true
			await this.updateData(dataPrevious)
			this.blockHistory = false
			this.redraw(true)
			console.log('<-- history_reverse:', JSON.stringify(this.history_reverse.map(state => state[0]['Test'])), this.history_reverse)
			console.log('<-- history:', JSON.stringify(this.history.map(state => state[0]['Test'])), this.history)
		}
	}

	// Redo undone action
	async redoChange() {
		if (this.history_reverse.length > 0) {
			const data = this.history_reverse.pop()
			this.history.push(data)
			this.blockHistory = true
			this.updateData(data)
			this.blockHistory = true
			this.redraw(true)
			console.log('--> history:', JSON.stringify(this.history.map(state => state[0]['Test'])))
			console.log('--> history_reverse:', JSON.stringify(this.history_reverse.map(state => state[0]['Test'])))
		}
	}

	// #endregion

	/////////////////////////////////
	// #region - Interaction

	onKeyDown(e) {
		const inputInFocus = e.target.tagName.toLowerCase() == 'input' || e.target.tagName.toLowerCase() == 'textarea'
		if (inputInFocus) return

		const $focusCell = document.querySelector('.tabulator-cell.focus')
		const $row = $focusCell ? $focusCell.closest('.tabulator-row') : null
		const row = $row ? this.getRow($row) : null

		// Copy focused cell content on cmd + C
		if ($focusCell && e.key == 'c' && (e.metaKey || e.ctrlKey)) {
			contextMenus.copyCell(null, null, $focusCell)
			return
		}

		if (e.key == 'Escape') {
			// Hide full cell content overlay
			const fullCellDisplayActive = this.hideFullCellContent()
			if (fullCellDisplayActive) return

			// Exit edit mode
			if (this.isEditMode()) {
				toggleEditMode(false, true)
				return
			}

			// Deselect rows
			if (this.getSelectedRows().length) {
				this.deselectRows(true)
				return
			}
		} else if ((e.metaKey || e.ctrlKey) && e.key == 'z') {
			// Undo / redo
			if (e.shiftKey) {
				this.redoChange()
			} else {
				this.undoChange()
			}
			return
		}

		// Move focus with arrow keys
		if ($focusCell) {
			if (e.key == 'ArrowLeft' || (e.shiftKey && e.key == 'Tab')) {
				this.moveFocus('left')
				e.preventDefault()
			} else if (e.key == 'ArrowRight' || e.key == 'Tab') {
				this.moveFocus('right')
				e.preventDefault()
			} else if (e.key == 'ArrowUp') {
				this.moveFocus('up')
				e.preventDefault()
			} else if (e.key == 'ArrowDown') {
				this.moveFocus('down')
				e.preventDefault()
			} else if (e.key == 'Escape') {
				$focusCell.classList.remove('focus')
				e.preventDefault()
			} else if (e.key == 'Enter') {
				if (this.isEditMode()) {
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
			} else if (this.isEditMode() && e.key.length == 1) {
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
	}

	// Homebrewed version of movableRows.
	// - - -
	// We had to build our own because the build-in movableRows
	// can't be made conditional, and it creates undesireable
	// side effects in edit mode. Specifically, it will interfere
	// with the textarea resize handle.
	// - - -
	// This comes at a price. We can't use the built-in download
	// function and we need to manually rearrange the rows according
	// to their UI order before we export data to CSV using our own
	// exporter. I think down the line, we should either rewrite this
	// as a module so we can plug into the built-in architecture for
	// rearranging rows, or we should submit a PR to Tabulator to make
	// their movableRows conditional.
	onCellMouseDown(clickEvent, cell) {
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
		function _onMouseUp(e) {
			document.removeEventListener('mousemove', _onMouseMove)
			document.removeEventListener('mouseup', _onMouseUp)

			// If no dragging occured, it's a click.
			if (!dragging) {
				table.onCellClick(e, cell)
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

	// The click evens is controlled from the onCellMouseDown() event.
	onCellClick(e, cell) {
		// Activate edit mode when cmd-clicking a cell.
		if (e.metaKey || e.ctrlKey) {
			toggleEditMode(true)
			return
		}

		// Display overlay with full content when cell is truncated.
		if (e.target.classList.contains('text-wrap')) {
			const $textWrap = e.target
			const $cell = $textWrap.closest('.tabulator-cell')
			const isTruncated = _isTruncated($textWrap)
			if (isTruncated) {
				this.displayFullCellContent($textWrap, $cell)
				// _expandCell($textWrap, $cell, cell)
				return
			}
		}

		// Set focus on clicked cell.
		const $focusCell = cell.getElement()
		this.setFocus($focusCell)

		// Select row
		this.onRowClick(e, cell.getRow())

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
			let $clone = $textWrap.cloneNode()
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
			$clone.remove() // From DOM
			$clone = null // From memory

			return isTruncated
		}

		// Copy styles from one element to another
		function _copyStyle($sourceElm, $targetElm) {
			const keys = ['width', 'font-size', 'line-height', 'padding']
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
	onRowClick(e, row) {
		// const currentRowIndex = row.getPosition()
		// - - -
		// We can't use getPosition because we're not using
		// the native movableRows. See table.onCellMouseDown().
		// So instead we have to figure out position from HTML:
		const $row = row.getElement()
		const currentRowIndex = Array.from($row.parentNode.querySelectorAll('.tabulator-row')).indexOf($row) + 1
		console.log(33, currentRowIndex)

		const lastSelectedRowIndex = this.lastSelectedRowIndex
		if (e.shiftKey && this.lastSelectedRowSelState != null) {
			// Select or deselect all rows between the last selected row and the current row
			const selectedRows = this.getSelectedRows()
			if (selectedRows.length) {
				let lowIndex = Math.min(lastSelectedRowIndex, currentRowIndex)
				let highIndex = Math.max(lastSelectedRowIndex, currentRowIndex)

				// When you select from bottom to top, we gotta include the highIndex
				// When you select from top to bottom, we gotta include the lowIndex
				if (lowIndex != lastSelectedRowIndex) {
					lowIndex -= 1
					highIndex -= 1
				}

				const toBeSelected = this.getRowsOrdered().slice(lowIndex, highIndex)
				if (this.lastSelectedRowSelState) {
					this.selectRow(toBeSelected)
					this.lastSelectedRowSelState = true
				} else {
					this.deselectRow(toBeSelected)
					this.lastSelectedRowSelState = false
				}
			}
		} else {
			// Toggle single row
			row.toggleSelect()
			this.lastSelectedRowSelState = this.getSelectedRows().includes(row)
		}
		this.lastSelectedRowIndex = currentRowIndex
		this.selectMode = this.getSelectedRows().length > 0
	}

	// Display overlay div that matches the cell's content and position,
	// but that expands to the bottom to display the full, untruncated content.
	displayFullCellContent($textWrap, $cell) {
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
	hideFullCellContent() {
		const $display = document.getElementById('display-full-text')
		if ($display) {
			$display.blur()
			return true
		}
		return false
	}

	// Deselect all rows, but ask user if it's more than 3
	deselectRows(soft) {
		const selectedRows = this.getSelectedRows()
		// Note: == true is on purpose, because we want to ignore the click event
		if (soft == true && selectedRows.length > 3) {
			if (confirm('Are you sure you want to deselect all rows?')) {
				this.deselectRow()
			}
		} else {
			this.deselectRow()
		}
	}

	// Set artificial cell focus
	setFocus($focusCell) {
		const $currentFocusCell = document.querySelector('.tabulator-cell.focus')
		if ($currentFocusCell) $currentFocusCell.classList.remove('focus')
		$focusCell.classList.add('focus')
	}

	// Unset artificial cell focus
	unsetFocus() {
		const $focusCell = document.querySelector('.tabulator-cell.focus')
		if ($focusCell) $focusCell.classList.remove('focus')
	}

	// Move artificial cell focus
	moveFocus(dir) {
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
	// #region - Tabulator overrides

	// Recalculate column widths when redrawing the table.
	redraw(bool) {
		// console.log('redraw', bool)
		if (!bool) {
			super.redraw()
			return
		}

		// Store current column widths.
		const colWidthsBefore = {}
		this.getColumns().forEach(col => {
			colWidthsBefore[col.getField()] = col.getWidth()
		})

		// Redraw
		super.redraw(bool)

		// Store new column widths.
		const colWidthsAfter = {}
		this.getColumns().forEach(col => {
			colWidthsAfter[col.getField()] = col.getWidth()
		})

		// console.log('colWidthsBefore:', colWidthsBefore)
		// console.log('colWidthsAfter:', colWidthsAfter)

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
			const isTampered = this.tamperedCols.includes(colName)
			// prettier-ignore
			const newWidth = isTampered
				? colWidthsBefore[colName]
				: Math.min(width, 400)
			this.getColumn(colName).setWidth(newWidth)
			// console.log(colName, isUntampered, '->', newWidth, width)
		})
	}

	// #endregion

	/////////////////////////////////
	// #region - Export

	// The native getRows() will return rows not in order because
	// we're not using the native movableRows. See table.onCellMouseDown().
	getRowsOrdered(selectedOnly) {
		let rows = this.getRows()

		// Rearrange data rows to match the order of the UI rows.
		const rowSelector = selectedOnly ? '.tabulator-row.tabulator-selected' : '.tabulator-row'
		const $rows = this.getRows()[0].getElement().closest('.tabulator-table').querySelectorAll(rowSelector)
		rows = Array.from($rows).map($row => {
			const $indexCol = $row.querySelector(`.tabulator-cell[tabulator-field=${this.options.index}]`)
			const rowIndex = +$indexCol.innerText - 1
			return rows[rowIndex]
		})
		return rows
	}

	// Prepare data for export.
	// - - -
	// This takes care of a few things that don't happen during UI manipulation:
	// 1. Remove deleted column fields from all rows
	// 2. Rearrange rows in the dataset to match the UI
	// - - -
	// Note: this is a little hacky – see onCellMouseDown for more info.
	getDataFinal(selectedOnly) {
		// let data = selectedOnly ? this.getSelectedData() : this.getData()
		let data = this.getData()
		const fields = this.getColumns().map(col => col.getField())

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
			const indexName = this.options.index
			data.forEach(row => {
				delete row[indexName]
			})
		}

		// Reorder data columns to match the order of the UI columns.
		data = data.map(row => {
			const newRow = {}
			fields.forEach(field => {
				newRow[field] = row[field]
			})
			return newRow
		})

		// Rearrange data rows to match the order of the UI rows.
		const rowSelector = selectedOnly ? '.tabulator-row.tabulator-selected' : '.tabulator-row'
		const $rows = this.getRows()[0].getElement().closest('.tabulator-table').querySelectorAll(rowSelector)
		data = Array.from($rows).map($row => {
			const $indexCol = $row.querySelector(`.tabulator-cell[tabulator-field=${this.options.index}]`)
			const rowIndex = +$indexCol.innerText - 1
			return data[rowIndex]
		})

		return data
	}

	// #endregion
}

Table.moduleName = 'custom'
