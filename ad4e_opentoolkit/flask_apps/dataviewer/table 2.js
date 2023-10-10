// This is a custom module that extends
// Tabulator with some extra functionality.
class Table extends Tabulator {
	constructor(element, options) {
		super(element, options)
		super.on('renderComplete', () => {
			this.storeCols()
		})
		super.on('columnResized', function () {
			this.element.classList.add('resized')
		})
	}

	// Store column default widths.
	storeCols() {
		this.defaultColumnWidths = {}
		this.getColumns().forEach(col => {
			this.defaultColumnWidths[col.getField()] = col.getWidth()
		})
	}

	// Reset column widths to default.
	resetCols() {
		this.getColumns().forEach(col => {
			col.setWidth(this.defaultColumnWidths[col.getField()])
			this.element.classList.remove('resized')
		})
	}

	toggleEditMode() {
		console.log(this.options.columns)
		// this.options.columns.forEach(col => {
		// 	col.editor = true
		// })
		// console.log(this.options.columns)
		// this.redraw(true)
	}
}

Table.moduleName = 'custom'
