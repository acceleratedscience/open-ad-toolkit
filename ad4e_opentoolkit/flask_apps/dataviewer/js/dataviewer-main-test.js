//define row context menu contents
var rowMenu = [
	{
		label: 'Test',
		action: function (e, row) {
			alert('test')
		},
	},
]

//define column header menu as column visibility toggle
var headerMenu = function () {
	return [
		{
			label: 'Test',
			action: function (e, row) {
				alert('test')
			},
		},
	]
}

// initialize table
var table = new Tabulator('#table', {
	rowContextMenu: rowMenu, //add context menu to rows
	data: [
		{ name: 'Foo', gender: 'Male' },
		{ name: 'Bar', gender: 'Female' },
	],
	columns: [
		{ title: 'Name', field: 'name', headerMenu: headerMenu },
		{ title: 'Gender', field: 'gender', headerMenu: headerMenu },
	],
})

const table1 = new Tabulator('#table-1', {
	rowContextMenu: rowMenu,
	data: [
		{ A: 'Foo', B: 'Bar' },
		{ A: 'Boo', B: 'Dar' },
	],
	columns: [
		{ title: 'A', field: 'A', headerMenu: headerMenu },
		{ title: 'B', field: 'B', headerMenu: headerMenu },
	],
})
