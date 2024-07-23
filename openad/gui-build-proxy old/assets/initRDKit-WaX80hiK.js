const o=(()=>{let i;return()=>(i||(i=new Promise((e,n)=>{window.initRDKitModule().then(t=>{window.RDKit=t,e(t)}).catch(t=>{n()})})),i)})();export{o as i};
