// Note: this file is organized in regions, which collapse in VS Code:
// - Collapse all regions:       hold cmd, then hit K followed by 1
// - Expand all regions:         hold cmd, then hit K followed by J
// - Collapse to any level:      first expand everything, then hold cmd, then hit K followed by any number

/**/
/**/
// #region - Main

class Table {
	constructor(data, headers) {
		if (!data) {
			console.error('No data provided to the TableControls class.')
			return
		}

		// Parse data
		this.data = JSON.parse(data)
		this.headers = JSON.parse(headers)
		this.colCount = this.data[0].length + 1 // +1 for checkbox column

		// Set column props
		const gridTemplateColumnsDefault = 'minmax(auto, 200px)'
		this.gridTemplateColumns = Array(this.colCount).fill(gridTemplateColumnsDefault)
		this.gridTemplateColumns[0] = '30px' // Checkbox column

		// Store DOM elements
		this.$table = document.getElementById('table')

		// Render HTML
		this.renderTable()

		// Initialize interaction
		this.initInteraction()
	}
}

// Render table HTML
Table.prototype.renderTable = function () {
	// Render header
	this.headers.unshift('') // Checkbox column
	this.headers.forEach((headerText, i) => {
		const $headerCell = document.createElement('div')
		$headerCell.classList.add('cell', 'header')
		if (i == 0) {
			$headerCell.classList.add('cb')
			headerText = '<input type="checkbox" id="checkAll">'
		}
		$headerCell.innerHTML = headerText
		if (i > 0) {
			$resizeHandle = document.createElement('div')
			$resizeHandle.classList.add('resize')
			$headerCell.append($resizeHandle)
		}
		this.$table.appendChild($headerCell)
	})

	// Render data
	this.data.forEach(row => {
		row.unshift('') // Checkbox column
		row.forEach((cellText, i) => {
			const $cell = document.createElement('div')
			$cell.classList.add('cell')
			if (i == 0) {
				$cell.classList.add('cb')
				cellText = '<input type="checkbox" name="checkRow">'
			} else {
				$cell.setAttribute('contenteditable', '')
			}
			$cell.innerHTML = cellText
			this.$table.appendChild($cell)
		})
	})

	// Set grid CSS
	this.setStyle()
}

// Insert or update a style tag with dynamic values.
Table.prototype.setStyle = function () {
	css = `
	#table {
		grid-template-columns: ${this.gridTemplateColumns.join(' ')};
	}
	#table > .cell:nth-child(${this.colCount}n) {
		border-right: solid 1px rgba(0, 0, 0, 0.1);
	}
	`
	let $style = document.getElementById('jsStyle')
	if (!$style) {
		$style = document.createElement('style')
		$style.setAttribute('id', 'jsStyle')
		document.head.append($style)
	}
	$style.innerHTML = css
}

// Initialize interaction
Table.prototype.initInteraction = function () {
	this.initResize()
	this.initSelect()
}

// #endregion

/**/
/**/
// #region - Resize columns

// Allow user to resize columns
Table.prototype.initResize = function () {
	// Bind handler to class scope.
	this.resizeColumn = this.resizeColumn.bind(this)

	// Start resizing - mouse down
	resizeHandles = document.querySelectorAll('#table .resize')
	resizeHandles.forEach(resizeHandle => {
		resizeHandle.addEventListener('mousedown', this.resizeStart.bind(this))
	})

	// Stop resizing - mouse up / window blur
	document.addEventListener('mouseup', this.resizeStop.bind(this))
	window.addEventListener('blur', this.resizeStop.bind(this))
}

// Start resizing.
Table.prototype.resizeStart = function (e) {
	// Check which column is being resized
	this.dragColNr = Array.from(e.target.parentNode.parentNode.children).indexOf(e.target.parentNode)

	// Store starting values
	this.startX = e.pageX
	this.startWidth = this.$table.querySelectorAll('.cell.header')[this.dragColNr].clientWidth

	// Fixate columns
	this.fixateColumns()

	// Prevent header text selection
	// this.dragging = true
	document.getElementById('table').classList.add('dragging')

	// Listen for dragging action
	document.addEventListener('mousemove', this.resizeColumn)
}

// Stop resizing
Table.prototype.resizeStop = function () {
	// Re-enable text selection
	// this.dragging = false
	document.getElementById('table').classList.remove('dragging')

	document.removeEventListener('mousemove', this.resizeColumn)
}

// Before resizing, we fixate all column widths
// to prevent the other columns from resizing.
Table.prototype.fixateColumns = function () {
	this.gridTemplateColumns.forEach((width, i) => {
		this.gridTemplateColumns[i] = `${this.$table.querySelectorAll('.cell.header')[i].clientWidth + 1}px`
	})
	this.setStyle()
}

// Resize column while dragging
Table.prototype.resizeColumn = function (e) {
	const dist = e.pageX - this.startX
	newWidth = Math.max(this.startWidth + dist, 30)
	this.gridTemplateColumns[this.dragColNr] = `${newWidth}px`
	this.$table.style.gridTemplateColumns = this.gridTemplateColumns.join(' ')
}

// #endregion

/**/
/**/
// #region - UI Options

Table.prototype.initSelect = function () {}

// Set selection
Table.prototype.selectRow = function (rowNr) {
	const $row = this.$table.querySelectorAll('.cell')[rowNr * this.colCount]
	$row.classList.toggle('selected')
}

// Set truncation
Table.prototype.setTruncate = function (trunc) {
	if (trunc == 'trunc-none') {
		this.$table.classList.remove('trunc-1', 'trunc-3')
	} else if (trunc == 'trunc-1') {
		this.$table.classList.remove('trunc-3')
		this.$table.classList.add('trunc-1')
	} else if (trunc == 'trunc-3') {
		this.$table.classList.remove('trunc-1')
		this.$table.classList.add('trunc-3')
	}
}

// #endregion
