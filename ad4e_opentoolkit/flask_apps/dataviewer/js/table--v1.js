// This is a custom module that extends
// Tabulator with some extra functionality.
class Table extends Tabulator {
	constructor(element, options) {
		// Modifications
		super(element, options)

		super.on('dataLoaded', () => {
			console.log('- - - dataLoaded - - -')
			if (this.addedIndexRow) return
			const result = this.addIndexCol(options.data)
			console.log('Result:', result)
			const { data, columns } = result
			console.log('@@@@', options.columns, [...columns])
			if (data && columns) {
				console.log('*', 1)
				this.data = data
				console.log('*', 2)
				this.columns = columns
				console.log('*', 3)
			}
		})

		super.on('tableBuilt', () => {
			console.log('*', 4)
			console.log(9999, !!this.data, !!this.columns, this.data, this.columns)
			if (this.data && this.columns) {
				console.log('> data:', this.data)
				console.log('> columns:', this.columns)
				this.setData(this.data)
				this.setColumns(this.columns)
				console.log(555, this.getColumns())
			}
		})
		// super.on('renderComplete', () => {
		// 	this.storeColWidths()
		// })
		// super.on('columnResized', () => {
		// 	this.element.classList.add('resized')
		// })

		// Variables
		this.editMode = false
		this.addedIndexCol = false
		this.colDefaultWidths = {}
		this.addedIndexRow = false

		// Parse custom options
		// this.init(options)
	}

	// Initialization
	init(options) {
		this.on('dataProcessed', () => {
			console.log('- - - dataProcessed - - -')
			// Turn on edit mode if hash is present
			if (options.editMode) {
				this.toggleEditMode(true)
			}
		})
	}

	// Trash
	// on(eventName, callback) {
	// 	console.log(333, eventName)
	// 	if (eventName != 'cellEditing') {
	// 		super.on(eventName, callback)
	// 	} else {
	// 		console.log(1234)
	// 	}
	// }

	// Store column default widths
	storeColWidths() {
		this.getColumns().forEach(col => {
			this.colDefaultWidths[col.getField()] = col.getWidth()
		})
	}

	// Add index column if not yet included
	addIndexCol(data) {
		console.log('A')
		if (this.addedIndexRow) return {}
		console.log('B')
		data = data ? data : this.getData()
		const rowSample = data[1]
		let hasIndex = !!rowSample['#']
		if (!hasIndex) {
			for (const key in rowSample) {
				if (key.toLowerCase() == 'index') {
					hasIndex = true
					break
				}
			}
		}
		if (!hasIndex) {
			// data.forEach((row, i) => {
			// 	row['#'] = i + 1
			// 	this.addedIndexRow = true
			// })

			const columns = this.getColumns()
			const colSample = columns[0]
			// columns[0] = {
			// 	...colSample,
			// 	sorter: 'number',
			// 	title: '#',
			// 	field: '#',
			// }
			console.log(11, data)
			console.log(22, columns)

			return { data, columns }
		} else {
			return {}
		}
	}

	// Reset column widths to default
	resetCols() {
		this.getColumns().forEach(col => {
			col.setWidth(this.colDefaultWidths[col.getField()])
			this.element.classList.remove('resized')
		})
	}

	// Toggle edit mode
	toggleEditMode(bool, revertChanges) {
		this.editMode = bool == undefined ? !this.editMode : bool
		if (this.editMode) {
			this.storeData()
			this.element.classList.add('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search + '#edit') // Add hash
		} else {
			if (revertChanges) this.revertData()
			this.element.classList.remove('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search) // Remove hash
		}
	}

	// Store copy of data so we can revert changes
	storeData() {
		this.dataBeforeEdit = this.getData()
		console.log('Store:', this.dataBeforeEdit)
	}

	// Revert changes on cancel
	revertData() {
		this.setData(this.dataBeforeEdit)
		console.log('Revert:', this.data, this.dataBeforeEdit)
		this.redraw(true)
	}
}

Table.moduleName = 'custom'
