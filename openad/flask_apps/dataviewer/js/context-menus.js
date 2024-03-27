// The content of the context menu depends on the state of the table.
// This class takes care of the logic and provides the menu items.

class ContextMenus {
	constructor() {
		this.header = this.contextMenuHeader.bind(this)
		this.cell = this.contextMenuCell.bind(this)
		this.separator = { separator: true }

		// Set on init
		this.table = null
		this.editAction = null
		this.saveEdits = null
		this.cancelEdits = null
	}

	init(table, { saveEdits, cancelEdits }) {
		this.table = table
		this.editAction = {
			label: 'Edit',
			action: this.editCell.bind(this),
			disabled: cell => {
				const table = cell.getTable()
				// When the dataset includes an index column, we set editable=False
				// for this column alone, after the initiation of table. Because of
				// that, we check for !isEditable OR isIndexColumn. See table.ensureUniqueIndex()
				const isIndexColumn = cell.getColumn().getField() == table.options.index
				const isEditable = cell.getColumn().getDefinition().editable
				if (isIndexColumn || !isEditable) return true
			},
		}
		this.saveEdits = saveEdits
		this.cancelEdits = cancelEdits
	}

	// Assembled context menu for the header.
	contextMenuHeader() {
		if (this.table.isEditMode()) return this.editModeActions()
		return [
			{
				label: 'Rename column',
				action: this.table.renameColumn,
			},
			{
				label: 'Delete column',
				action: this.table.deleteColumn,
			},
		]
	}

	// Assembled context menu for a cell.
	contextMenuCell(e, cell) {
		if (this.table.isEditMode()) return this.editModeActions()
		let menuItems = []

		if (cell.getRow().isSelected()) {
			// Display selection actions first, if row is selected.
			// prettier-ignore
			menuItems = [
                ...this.selectedActions(),
                this.separator,
                this.editAction,
                { label: 'Row', menu: this.rowActions() },
                { label: 'Cell', menu: this.cellActions() }
            ]
		} else {
			// Default menu
			// prettier-ignore
			menuItems = [
                this.editAction,
                this.separator,
                ...this.rowActions(),
                this.separator,
                ...this.cellActions()
            ]
		}

		return menuItems
	}

	/**
	 * Context menu items organized by category.
	 */

	cellActions() {
		return [
			{
				label: 'Copy cell',
				action: this.copyCell,
			},
		]
	}

	rowActions() {
		return [
			{
				label: 'Delete row',
				action: (e, cell) => {
					const row = cell.getRow()
					row.delete()
				},
			},
			{
				label: 'Copy row',
				action: this.copyRow,
			},
		]
	}

	selectedActions() {
		return [
			{
				label: 'Delete selected',
				action: this.deleteSelected.bind(this),
			},
			{
				label: 'Keep selected',
				action: this.keepSelected.bind(this),
			},
			{
				label: 'Copy selected',
				action: () => this.copyData.bind(this)(true),
			},
			{
				label: 'Download selected',
				action: () => this.downloadData.bind(this)(true),
			},
		]
	}

	editModeActions() {
		return [
			{
				label: 'Save',
				action: this.saveEdits,
			},
			{
				label: 'Cancel',
				action: this.cancelEdits,
			},
		]
	}

	/**
	 * Cell actions
	 */

	// Toggle edit mode and select cell.
	editCell(e, cell) {
		const field = cell.getField()
		if (field == this.table.options.index) {
			let message = `You can't edit the '${field}' column as it's used for identifying the rows.` // Note: repeat
			// To do: we should provide an option to create an additional index column.
			// if (!this.table.addedIndex) {
			// 	message += " In order for it to be editable, you can change the column name to something other than 'index'."
			// }
			alert(message)
		} else {
			toggleEditMode(true)
			const $row = cell.getElement()
			$row.click()
		}
	}

	copyCell(e, cell, $cell) {
		$cell = $cell || cell.getElement()
		navigator.clipboard.writeText($cell.innerText)

		// Blink
		$cell.classList.add('copied')
		setTimeout(() => {
			$cell.classList.remove('copied')
		}, 600)
	}

	/**
	 * Row actions
	 */

	copyRow(e, cell) {
		const row = cell.getRow()
		const rowData = row.getData()
		const text = Object.values(rowData).join('\t')
		navigator.clipboard.writeText(text)

		// Blink
		const $row = row.getElement()
		$row.classList.add('copied')
		setTimeout(() => {
			$row.classList.remove('copied')
		}, 600)
	}

	/**
	 * Selection actions
	 */

	deleteSelected() {
		const selectedIndices = this.table.getSelectedRows().map(row => row.getIndex())
		this.table.getRows().forEach(row => {
			if (selectedIndices.includes(row.getIndex())) {
				row.delete()
			}
		})
		this.table.deselectRows()
	}

	keepSelected() {
		const selectedIndices = this.table.getSelectedRows().map(row => row.getIndex())
		this.table.getRows().forEach(row => {
			if (!selectedIndices.includes(row.getIndex())) {
				row.delete()
			}
		})
		this.table.deselectRows()
	}

	// Copy selected rows buy default, or all rows if all is true.
	copyData(selectedOnly) {
		const data = this.table.getDataFinal(selectedOnly)
		const fields = this.table.getColumns().map(col => col.getField())
		const textHeader = Object.values(fields).join('\t')
		const textRows = data.map(row => Object.values(row).join('\t')).join('\n') // prettier-ignore
		const text = textHeader + '\n' + textRows
		navigator.clipboard.writeText(text)

		// Blink
		const rows = selectedOnly ? this.table.getSelectedRows() : this.table.getRows()
		rows.forEach(row => {
			const $row = row.getElement()
			$row.classList.add('copied')
			setTimeout(() => {
				$row.classList.remove('copied')
			}, 600)
		})
	}

	// Download data, everything or selected only.
	downloadData(selectedOnly) {
		const data = this.table.getDataFinal(selectedOnly)

		const filename = prompt('Filename (.csv)', 'mydata')
		if (!filename) return

		// Compile text.
		const text = _compileCSVText.bind(this)()

		// Create a blob from the text.
		const blob = new Blob([text], { type: 'text/plain' })

		// Create anchor element.
		const a = document.createElement('a')
		a.href = URL.createObjectURL(blob)
		a.download = filename + '.csv'

		// Append the anchor to the document, trigger
		// a click (to download), then remove it.
		document.body.appendChild(a)
		a.click()
		document.body.removeChild(a)

		//
		//

		function _compileCSVText() {
			const headers = this.table.getColumns().map(col => col.getField())
			const textHeader = Object.values(headers).join(',')
			const textRows = data.map(row => Object.values(row).join(',')).join('\n') // prettier-ignore
			const text = textHeader + '\n' + textRows
			return text
		}
	}
}
