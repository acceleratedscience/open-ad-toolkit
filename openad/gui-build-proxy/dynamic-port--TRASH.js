// This file replaces the base URL of our app with one
// that includes the correct port number, which is dynamic.

// const port = process.env.PORT || '8024';
const base = `${location.pathname}proxy/`;
__webpack_public_path__ = base;