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

// Create table
const table = new Table('#table', {
	data,
	columns,
	resizableRows: true,
	// rowHeight: 40,
	reactiveData: true,
	editMode: window.location.hash == '#edit',
})

// Parse data & create columns.
function parseColumns(data) {
	const dataTypes = {}
	const keys = {}

	// Store the data types of each column of each row.
	data.forEach((row, i) => {
		Object.entries(row).map(([key, val]) => {
			val = val ? val.toString() : ''
			if (!dataTypes[key]) dataTypes[key] = []
			if (isDate(val)) {
				dataTypes[key].push('date')
				row[key] = moment(val).format('YYYY-MM-DD')
			} else if (val.match(/^http(s)?:\/\//)) {
				dataTypes[key].push('url')
				row[key] = val.replace(/^http(s)?:\/\/([a-zA-Z0-9$-_.+!*'(),/&?=:%]+?)(\/)?$/, '$2')
			} else {
				dataTypes[key].push(typeof val)
			}
		})
	})

	// Compare the data types for each row, and if
	// they all match, set the column type to that.
	Object.entries(dataTypes).map(([key, types]) => {
		if (types.every((type, i, arr) => type == arr[0])) {
			keys[key] = types[0]
		} else {
			keys[key] = 'string'
		}
	})

	// Create columns
	const columns = Object.entries(keys).map(([key, type]) => {
		col = {
			sorter: type,
			title: key,
			field: key,
			editor: true,
			editable: isEditMode,
			maxWidth: 500,
		}
		if (type == 'date') {
			col.sorterParams = {
				format: 'yyyy-MM-dd',
				alignEmptyValues: 'top',
			}
		} else if (type == 'number') {
			col.sorterParams = {
				thousandSeparator: ',',
				decimalSeparator: '.',
				alignEmptyValues: 'top',
			}
		} else if (type == 'url') {
			col.sorter = 'string'
			col.formatter = 'link'
			col.formatterParams = { target: '_blank' }
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

	console.log(keys)
	console.log(data)
	console.log(columns)
	return columns
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
