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

let data = [
	{ id: 1, name: 'Billy Bob', progress: '12', gender: 'male', height: 1, col: 'red', dob: '30/11/1984', driver: 1, color: 'red' },
	{ id: 2, name: 'Mary May', progress: '1', gender: 'female', height: 2, col: 'blue', dob: '14/05/1982', driver: true, color: 'green' },
	{ id: 3, name: 'Christine Lobowski', progress: '42', height: 0, col: 'green', dob: '22/05/1982', driver: 'true', color: 'yellow' },
	{ id: 4, name: 'Brendon Philips', progress: '125', gender: 'male', height: 1, col: 'orange', dob: '01/08/1980', color: 'blue' },
	{ id: 5, name: 'Margret Marmajuke', progress: '16', gender: 'female', height: 5, col: 'yellow', dob: '31/01/1999', color: 'orange' },
]

const table = new Table('#table', {
	data: data,
	reactiveData: true,
	editMode: window.location.hash == '#edit',
	columns: [
		{ editor: true, editable: isEditMode, title: 'Name', field: 'name', sorter: 'string', width: 200 },
		{ editor: true, editable: isEditMode, title: 'Progress', field: 'progress', sorter: 'number', formatter: 'progress' },
		{ editor: true, editable: isEditMode, title: 'Gender', field: 'gender', sorter: 'string' },
		{ editor: true, editable: isEditMode, title: 'Rating', field: 'rating', formatter: 'star', hozAlign: 'center', width: 100 },
		{ editor: true, editable: isEditMode, title: 'Favourite Color', field: 'col', sorter: 'string' },
		{ editor: true, editable: isEditMode, title: 'Date Of Birth', field: 'dob', sorter: 'date', hozAlign: 'center' },
		{ editor: true, editable: isEditMode, title: 'Driver', field: 'car', hozAlign: 'center', formatter: 'tickCross', sorter: 'boolean' },
		{ editor: true, editable: isEditMode, title: 'Color', field: 'color', hozAlign: 'right', sorter: 'string' },
	],
})

function isEditMode(cell, a) {
	return table.editMode
}
