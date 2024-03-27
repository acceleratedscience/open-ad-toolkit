render3dMol()
document.addEventListener('DOMContentLoaded', function () {
	// Enrich molecule data with PubChem.
	if (!document.getElementById('data-cid')) {
		enrichMolData()
	}

	// Set height for synonyms & parameter blocks.
	setHeightSynonyms()
	setHeightParameters()

	// Enable expandable sections
	enableToggles()
})

//
//

/**
 * Connect to the /enrich API endpoint to fetch
 * additional data about the molecule.
 */
async function enrichMolData() {
	try {
		// Loading title
		document.getElementById('data-name').classList.add('loading')

		// Display fetching message
		document.getElementById('fetching-pubchem').innerHTML = document.getElementById('fetching-pubchem').getAttribute('data-text')
		document.getElementById('fetching-pubchem').classList.remove('error')

		// Fetch
		const response = await fetch(`/enrich`, {
			method: 'POST',
			body: document.getElementById('data-inchi').innerText,
		})

		if (response.ok) {
			// Update HTML
			const html = await response.text()
			if (html) {
				const dump = document.createElement('html')
				dump.innerHTML = html
				grid = dump.querySelector('#grid')
				document.getElementById('grid').innerHTML = grid.innerHTML
			}
		} else {
			// Handle errors
			const err = 'Failed to connect'
			document.getElementById('fetching-pubchem').innerHTML = err + ' - <a href="#" onclick="enrichMolData()">retry</a>'
			document.getElementById('fetching-pubchem').classList.add('error')
		}

		// Remove loading from molecule name.
		document.getElementById('data-name').classList.remove('loading')

		// Re-initialize UI elements.
		enableToggles()
		setHeightSynonyms()
		setHeightParameters()
	} catch (err) {
		console.error('Something went wrong fetching the enriched molecule data.', err)
	}
}

// Enable expandable sections
function enableToggles() {
	document.querySelectorAll('.toggle-expand').forEach($elm => {
		$elm.removeEventListener('click', _toggleExpand)
		$elm.addEventListener('click', _toggleExpand)
	})
}

function _toggleExpand(e) {
	e.preventDefault()
	e.target.classList.toggle('expand')
}

// Set height for synonyms wrapper.
// This makes sure they are divided in 4 neat columns.
function setHeightSynonyms() {
	const $wrap = document.querySelector('#synonyms .synonyms-wrap')
	const count = $wrap.querySelectorAll('div').length
	const height = Math.ceil(count / 4) * 22
	$wrap.setAttribute('style', `max-height: ${height}px`)

	// Set cloak min height - TO DO FIX THIS
	// const $cloak = document.querySelector('#synonyms .cloak')
	// const minHeight = Math.min(110, height)
	// $cloak.setAttribute('style', `min-height: ${minHeight}px`)

	// Display number of synonyms in the toggle.
	const $toggle = document.querySelector('#synonyms .toggle-expand')
	if (count) {
		$toggle.classList.remove('hide')
		$toggle.querySelector('span').innerText = ' ' + count
	}
}

// Set height for parameters wrapper.
// This makes sure they are divided in 4 neat columns.
function setHeightParameters() {
	const $wrap = document.querySelector('#parameters .param-wrap')
	const count = $wrap.children.length
	const height = Math.ceil(count / 3) * 22
	$wrap.setAttribute('style', `max-height: ${height}px`)
}

/**
 * Render the molecule in 3D.
 */
function render3dMol() {
	const $wrap = document.querySelector('#mol-render .mol-3d')
	const inchi = document.getElementById('data-inchi').innerText
	// _render3dMol_Jmol($wrap, inchi)
	_render3dMol_3DMol_UNUSED($wrap)
}

// Using the 3DMol library - https://3dmol.org
// Not great because it renders 3D structures flat.
function _render3dMol_3DMol_UNUSED($wrap) {
	// let config = { backgroundColor: '#333' }
	const config = { backgroundColor: '#000' }
	const viewer = $3Dmol.createViewer($wrap, config)

	const mol_sdf = document.getElementById('mol-data').getAttribute('data-sdf')

	viewer.addModel(mol_sdf, 'sdf')

	// Check styling options: https://3dmol.org/tests/example.html
	// More advances visualizing: https://3dmol.org/viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick&select=chain:B&style=line;sphere&select=resi:19,23,26;chain:B&style=stick;sphere&select=resi:19,23,26;chain:B&labelres=backgroundOpacity:0.8;fontSize:14
	viewer.setStyle({}, { stick: {} })
	// viewer.setStyle({}, { sphere: {} })
	// viewer.setStyle({ hetflag: false }, { cartoon: {} })
	viewer.addSurface($3Dmol.SurfaceType.MS, { map: { prop: 'partialCharge', scheme: new $3Dmol.Gradient.RWB(-0.6, 0.6) }, opacity: 0.85 }, { chain: 'B' }, { chain: 'B' })
	viewer.addSurface($3Dmol.SurfaceType.VDW, {}, { hetflag: false, chain: 'A' }, { hetflag: false, chain: 'A' })
	viewer.zoomTo()
	viewer.render()
	_colorSS(viewer)

	function _colorSS(viewer) {
		var m = viewer.getModel()
		m.setColorByFunction({}, function (atom) {
			// console.log(atom, atom.elem, atom.capDrawn, atom.hetflag, atom.bonds, atom.bondOrder)
			if (atom.elem == 'O') return 'red'
			else if (atom.elem == 'C') return 'yellow'
			else return 'gray'
		})
		viewer.render()
	}
}

// Using the Jmol library - https://jmol.sourceforge.net
// Better but still not great because it relies on 30MB library files.
function _render3dMol_Jmol($wrap, inchi) {
	const molInfo = {
		width: '100%',
		height: '100%',
		debug: false,
		j2sPath: 'app/jsmol/j2s',
		// js2Path: 'https://cdn.jsdelivr.net/npm/@ruiming0114/jsmol@latest/j2s', // Doesn't work because the folder is rendered as HTML
		color: '0xC0C0C0',
		disableJ2SLoadMonitor: true,
		disableInitialConsole: true,
		appletLoadingImage: 'none',
		addSelectionOptions: false,
		serverURL: 'app/jsmol/jsmol.php',
		use: 'HTML5',
		readyFunction: null,
		// script: 'load $' + smiles,
		script: 'load "$' + inchi + '"; background "#f6f6f6";', //  zoom 70;
	}

	$wrap.innerHTML = Jmol.getAppletHtml('jmolApplet1', molInfo)
}

function _parsePubChem(pubchem_data) {
	const params = {}

	const sections = pubchem_data.Record.Section
	sections.forEach(section => {
		if (section.TOCHeading == 'Names and Identifiers') {
			const subSections = section.Section
			subSections.forEach(subSection => {
				if (subSection.TOCHeading == 'Molecular Formula') {
					params[subSection.TOCHeading] = subSection.Information[0].Value.StringWithMarkup[0].String
				} else if (subSection.TOCHeading == 'Synonyms') {
					const subSubSections = subSection.Section
					subSubSections.forEach(subSubSection => {
						if (subSubSection.TOCHeading == 'Depositor-Supplied Synonyms') {
							params[subSection.TOCHeading] = subSubSection.Information[0].Value.StringWithMarkup.map(obj => obj.String)
						}
					})
				}
			})
		} else if (section.TOCHeading == 'Chemical and Physical Properties') {
			const subSections = section.Section
			subSections.forEach(subSection => {
				const propertyType = subSection.TOCHeading // 'Computed Properties' / 'Experimental Properties' / 'SpringerMaterials Properties' / etc --> will become table column
				if (subSection.TOCHeading == 'SpringerMaterials Properties') {
					// Ignore, this contains links, not values.
				} else {
					const subSubSections = subSection.Section
					if (subSubSections) {
						subSubSections.forEach(subSubSection => {
							const key = subSubSection.TOCHeading
							let value = ''
							if ('StringWithMarkup' in subSubSection.Information[0].Value) {
								value = subSubSection.Information[0].Value.StringWithMarkup[0].String
							} else if ('Number' in subSubSection.Information[0].Value) {
								value = subSubSection.Information[0].Value.Number[0]
							}
							if ('Unit' in subSubSection.Information[0].Value) {
								value += ' ' + subSubSection.Information[0].Value.Unit
							}
							params[key] = value
						})
					} else {
						console.error('Dev Warning: Subsection is not following usual structure, needs its own parser:', subSection.TOCHeading)
					}
				}
			})
		} else if (section.TOCHeading == 'Spectral Information') {
			// Todo
		} else if (section.TOCHeading == 'Related Records') {
			// Todo
		} else if (section.TOCHeading == 'Chemical Vendors') {
			// Ignore
		} else if (section.TOCHeading == 'Safety and Hazards') {
			// Todo
		} else if (section.TOCHeading == 'Literature') {
			// External table
		} else if (section.TOCHeading == 'Patents') {
			// Todo
		} else if (section.TOCHeading == 'Biological Test Results') {
			// External table
		} else if (section.TOCHeading == 'Taxonomy') {
			// External table
		} else if (section.TOCHeading == 'Classification') {
			// Todo
		}
	})
	return params
}
