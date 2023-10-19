document.getElementById('btn-save').addEventListener('click', submitData)

// Create a new XMLHttpRequest object
var xhr = new XMLHttpRequest()

function submitData() {
	// Define the request method and URL.
	xhr.open('POST', '/submit', true)

	// Set up a callback function to handle the response.
	xhr.onload = function () {
		if (xhr.status === 200) {
			// Success
			console.log('success')
			// alert(`You submitted: ${xhr.response}`)
			// window.location.href = `/success?someInput=${xhr.response}`
		} else {
			// Error
			console.log('error')
			// alert('Submit request failed with status code ' + xhr.status)
		}
	}

	// Send the request
	xhr.send(JSON.stringify(data))
}
