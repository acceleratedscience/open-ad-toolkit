render3dMol()
populateMolData(document.getElementById('data-inchi').innerText)

function render3dMol() {
	const $wrap = document.querySelector('#mol-render .mol-3d')
	// let config = { backgroundColor: '#333' }
	const config = { backgroundColor: '#000' }
	const viewer = $3Dmol.createViewer($wrap, config)

	const mol_sdf = document.getElementById('mol-data').getAttribute('data-sdf')

	viewer.addModel(mol_sdf, 'sdf')

	// Check styling options: https://3dmol.org/tests/example.html
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

// Update the DOM with PubChem data.
async function populateMolData(inchi) {
	let mol_json = document.getElementById('mol-data').getAttribute('data-json')
	mol_json = JSON.parse(mol_json)
	const cid = await fetchMolCID(inchi)
	if (cid) {
		mol_json['cid'] = cid
		const pubchem_data = await fetchMolData(cid)
		const pubchem_data_simple = parsePubChem(pubchem_data)
		mol_json.pubChem = pubchem_data_simple

		/**
		 * Update the DOM
		 */

		// Molecule name
		const currentName = document.getElementById('data-name').getAttribute('data-val')
		if (!currentName && 'RecordTitle' in pubchem_data.Record) {
			document.getElementById('data-name').innerText = pubchem_data.Record.RecordTitle
		}

		// CID
		document.getElementById('data-cid').innerText = cid
		document.getElementById('data-cid').setAttribute('href', 'https://pubchem.ncbi.nlm.nih.gov/compound/' + cid)

		// Parameters
		Object.entries(pubchem_data_simple).forEach(([key, val]) => {
			if (Array.isArray(val)) {
				const showMoreLink = val.length > 10 ? ` <a href="#" class="show-more"></a>` : ''
				document.getElementById('parameters').insertAdjacentHTML('beforeend', `<div><b>${key}: </b><span class="array"><span>${val.join('</span><span>, ')}</span>${showMoreLink}</span></div>`)
			} else {
				document.getElementById('parameters').insertAdjacentHTML('beforeend', `<div><b>${key}: </b>${val}</div>`)
			}
		})

		// Show-more links (masking long arrays)
		// To see this in action, check ibuprofen - show mol 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'
		document.querySelectorAll('.show-more').forEach($elm => {
			$elm.addEventListener('click', e => {
				e.preventDefault()
				$elm.closest('.array').classList.toggle('expand')
			})
		})

		// JSON data
		document.querySelectorAll('pre')[0].innerHTML = JSON.stringify(mol_json, null, '\t')
		document.querySelectorAll('pre')[1].innerHTML = JSON.stringify(pubchem_data, null, '\t')
	}
}

function parsePubChem(pubchem_data) {
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

// Fetch the CID from PubChem.
async function fetchMolCID(inchi) {
	try {
		// https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial#section=Special-Characters-in-the-URL
		const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=${encodeURIComponent(inchi)}`)
		if (!response.ok) {
			throw new Error('fetchMolCID() PubChem PUG REST API response is not ok')
		}
		const data = await response.json()
		const cid = 'IdentifierList' in data && 'CID' in data.IdentifierList ? data.IdentifierList.CID[0] : null
		if (!cid) {
			const err = Fault in data && Details in data.Fault && output in data.Fault.Details ? '\n' + data.Fault.Details.output : ''
			throw new Error('fetchMolCID() failed to extract cid from API response', err)
		}
		return cid
	} catch (err) {
		console.error(err)
	}
}

// Fetch all molecule data from PubChem.
async function fetchMolData(cid) {
	try {
		// https://pubchem.ncbi.nlm.nih.gov/docs/pug-view
		const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/${cid}/JSON`)
		if (!response.ok) {
			throw new Error('fetchMolData() PubChem PUG View API response is not ok')
		}
		const data = await response.json()
		return data
	} catch (err) {
		console.error(err)
	}
}
