import{b as m,z as M,M as S,N as h,e as a,A as o,O as I,P as u,K as v,L as w,d as P,h as d,i as f,t as p,G as b}from"./index-Bu5CNeR8.js";const y=m(),A=M(),C=100,E=["name","isomeric_smiles","molecular_formula"],F=["molecular_weight"];function g(){return{_context:null,_disableUpdate:!0,_hasChanges:!1,_cacheId:null,_mols:null,_total:0,_resultCount:0,_searchStr:"",_searchMode:"text",_page:1,_pageSize:C,_sort:"",_highlight:"",_focus:null,_sel:[],_allIndices:[],_matchingIndices:[],_availIdentifiers:[],_showIdentifiers:E,_availProps:[],_showProps:F}}const B=S("molGridStore",{state:()=>g(),getters:{context(){return this._context},active(){return!!this._cacheId},hasChanges(){return this._hasChanges},cacheId(){return this._cacheId},mols(){return this._mols},molSmiles(){return this._mols.map(e=>e.identifiers.isomeric_smiles||e.identifiers.canonical_smiles||e.identifiers.smiles)},resultCount(){return this._resultCount},total(){return this._total},searchStr(){return this._searchStr},searchMode(){return this._searchMode},page(){return this._page},pageTotal(){return Math.ceil(this._resultCount/this._pageSize)},pageSize(){return this._pageSize},sort(){return this._sort},highlight(){return this._highlight},focus(){return this._focus},sel(){return this._sel},hasSel(){return this._sel.length>0},allIndices(){return this._allIndices},matchingIndices(){return this._matchingIndices},availIdentifiers(){return this._availIdentifiers},showIdentifiers(){return this._showIdentifiers},availProps(){return this._availProps},showProps(){return this._showProps.filter(e=>this.availProps.includes(e))}},actions:{setContext(e){this._context=e},setHasChanges(e){console.log(">> ",e),this._hasChanges=e},async updateMols(e=!1){const s=e?{...h.currentRoute.value.query}:this._setUrlQuery();a(o.getMolset(this._cacheId,s),{onSuccess:t=>{this.setMolset(t)},onError:t=>{console.log("Error in getMolset()",t)}})},_setUrlQuery(){const e={...h.currentRoute.value.query};this.searchStr?e.search=this.searchStr:delete e.search,this.searchMode=="smarts"?e.smarts="1":delete e.smarts,this.page==1?delete e.page:e.page=String(this.page),this.sort?e.sort=this.sort:delete e.sort;const s=I(e),t=h.currentRoute.value.path+s;return h.push(t),e},async parseUrlQuery(){const e=h.currentRoute.value.query;this._disableUpdate=!0;let s=!1;e.search&&this._searchStr!=String(e.search)&&(this._searchStr=e.search||"",s=!0),this._searchMode!=(e.smarts?"smarts":"text")&&(this._searchMode=e.smarts?"smarts":"text",s=!0),e.page&&this._page!=Number(e.page)&&(this._page=Number(e.page),s=!0),this._sort!=(e.sort||"")&&(this._sort=e.sort||"",s=!0),await u(),this._disableUpdate=!1,s&&this.updateMols(!0)},async setMolset(e){this._disableUpdate=!0,this._cacheId=e.cacheId,this._mols=e.mols,this._total=e.total,this._resultCount=e.resultCount,this._searchStr=e.searchStr,this._searchMode=e.searchMode,this._sort=e.sort??"",this._allIndices=e.allIndices,this._matchingIndices=e.matchingIndices,this._pageSize=e.pageSize,this._page=e.page;const s=e.mols[0];s.properties.isomeric_smiles||(s.properties.canonical_smiles?(this.setIdentifier("canonical_smiles",!0),this.setIdentifier("isomeric_smiles",!1)):s.properties.smiles&&(this.setIdentifier("smiles",!0),this.setIdentifier("isomeric_smiles",!1))),e.searchMode=="smarts"&&(this._highlight=e.searchStr),e.mols.forEach(t=>{Object.keys(t.properties).forEach(i=>{!this._availProps.includes(i)&&t.properties[i]!=null&&this._availProps.push(i)}),Object.keys(t.identifiers).forEach(i=>{!this._availIdentifiers.includes(i)&&t.identifiers[i]!=null&&this._availIdentifiers.push(i)})}),await u(),this._disableUpdate=!1},removeMols(e){a(o.removeFromMolset(this._cacheId,e,h.currentRoute.value.query),{onSuccess:s=>{this.setMolset(s),this.deselectAll(),this.setHasChanges(!0)},onError:s=>{console.log("Error in removeFromMolset()",s)}})},async keepMols(e){const s=this._matchingIndices.filter(t=>!e.includes(t));this.removeMols(s)},replaceMolInMolset(e,s,t){return new Promise((i,r)=>{a(o.replaceMolInMolset(e,s,t,this._cacheId),{onSuccess:()=>{this.setHasChanges(!0),i(!0)},onError:n=>r(n)})})},updateMolset(){return new Promise((e,s)=>{a(o.updateMolset(y.path,this._cacheId),{onSuccess:()=>{this.setHasChanges(!1),e(!0)},onError:t=>{console.log("Error in updateMolset()",t),s(!0)}})})},_maybeShowInvalidMolsModal(e,{callback:s,destinationPath:t,options:i}){if(!e||typeof e!="object")return;const{data:r}=e;if(r.invalidMols){console.log(`The following ${r.invalidMols.length} molecules are invalid:`);const n=r.invalidMols.map(l=>{if(l){const c=l.identifiers.name?l.identifiers.name+"<br>":"",_=l.identifiers.isomeric_smiles||l.identifiers.canonical_smiles||l.identifiers.smiles;return console.log(l.index,_),`<li><b>${l.index}</b>: ${c}${_}</li>`}else return null}).join(`
`);return new Promise(l=>{A.alert(`The following molecule${r.invalidMols.length>1?"s":""} could not be processed by RDKit and will be removed:</p><ul>${n}</ul>`,{title:"Invalid molecules detected",html:!0,size:"md",primaryBtn:"Continue",secondaryBtn:!0,onSubmit:()=>{const c={...i};c.removeInvalidMols=!0,s(t,c),l(!0)}})})}},updateMolset_mymols(){return new Promise((e,s)=>{a(o.updateMolset_mymols(this._cacheId),{onSuccess:()=>e(!0),onError:()=>s(!0)})})},updateMolset_result(){return new Promise((e,s)=>{a(v.updateResult_molset(this._cacheId),{onSuccess:()=>e(!0),onError:()=>s(!0)})})},updateMolset_dataframe(e){return new Promise((s,t)=>{a(w.updateDataframe_molset(e,this._cacheId),{onSuccess:()=>s(!0),onError:()=>t(!0)})})},saveMolsetAsJSON(e,{newFile:s=!1,force:t=!1}={}){return new Promise((i,r)=>{a(o.saveMolsetAsJSON(e,this._cacheId,s,t),{onSuccess:()=>i(!0),onError:()=>r(!0)})})},async saveMolsetAsSDF(e,s){try{return await this._saveMolsetAsSDF(e,s)}catch(t){return this._maybeShowInvalidMolsModal(t,{callback:this.saveMolsetAsSDF,destinationPath:e,options:s}),Promise.resolve(!1)}},_saveMolsetAsSDF(e,{removeInvalidMols:s=!1,newFile:t=!1,force:i=!1}={}){return new Promise((r,n)=>{a(o.saveMolsetAsSDF(e,this._cacheId,s,t,i),{onSuccess:()=>r(!0),onError:l=>n(l)})})},async saveMolsetAsCSV(e,s={}){try{return await this._saveMolsetAsCSV(e,s)}catch(t){return this._maybeShowInvalidMolsModal(t,{callback:this.saveMolsetAsCSV,destinationPath:e,options:s}),Promise.resolve(!1)}},_saveMolsetAsCSV(e,{newFile:s=!1,force:t=!1}={}){return new Promise((i,r)=>{a(o.saveMolsetAsCSV(e,this._cacheId,s,t),{onSuccess:()=>i(!0),onError:()=>r(!0)})})},async saveMolsetAsSmiles(e,s={}){try{return await this._saveMolsetAsSmiles(e,s)}catch(t){return this._maybeShowInvalidMolsModal(t,{callback:this.saveMolsetAsSmiles,destinationPath:e,options:s}),Promise.resolve(!1)}},_saveMolsetAsSmiles(e,{newFile:s=!1,force:t=!1}){return new Promise((i,r)=>{a(o.saveMolsetAsSmiles(e,this._cacheId,s,t),{onSuccess:()=>i(!0),onError:()=>r(!0)})})},setSearchStr(e){this._searchStr=e||"",this.searchMode=="smarts"&&this.setHighlight(e),this._disableUpdate||this.updateMols()},setSearchMode(e){e=="text"?this.setHighlight(""):e=="smarts"&&this.setHighlight(this._searchStr),this._searchMode=e,console.log("setSearMode %%",e),this._disableUpdate||this.updateMols()},setPage(e){console.log("setPage"),this._page=e,this._disableUpdate||this.updateMols()},setPageSize(e){this._pageSize=e,this._disableUpdate||(this.updateMols(),this.resetPagination())},resetPagination(){this._page=1},setSort(e){console.log("setSort",e),e?this._sort=e:this._sort="",this._disableUpdate||this.updateMols()},getDisplayIndex(e){return this._matchingIndices.indexOf(e)+1},setHighlight(e){this._highlight=e},setFocus(e){this._focus=e},unsetFocus(){this._focus=null},setSel(e){this._sel=e},addSel(e){this._sel=[...new Set([...this._sel,...e])]},removeSel(e){this._sel=this._sel.filter(s=>!e.includes(s))},toggleSel(e){this._sel.includes(e)?this._sel=this._sel.filter(s=>s!==e):this._sel.push(e)},deselectAll(e=!1){e&&this._sel.length>3?confirm(`Deselect ${this._sel.length} molecules?`)&&(this._sel=[]):this._sel=[]},toggleIdentifier(e){this._showIdentifiers.includes(e)?this._showIdentifiers=this._showIdentifiers.filter(s=>s!==e):this._showIdentifiers.push(e)},setIdentifier(e,s){s===!1?this._showIdentifiers=this._showIdentifiers.filter(t=>t!==e):s===!0&&this._showIdentifiers.push(e)},toggleProp(e){e=e.replace(/^-/,""),e!="name"&&(this._showProps.includes(e)?this._showProps=this._showProps.filter(s=>s!==e):this._showProps.push(e))},enableProp(e){e=e.replace(/^-/,""),e!="name"&&(this._showProps.includes(e)||this._showProps.push(e))},async clear(){this._disableUpdate=!0,this.clearWorkingCopy(),Object.assign(this,g()),await u(),this._disableUpdate=!1},clearWorkingCopy(){a(o.clearMolsetWorkingCopy(this._cacheId),{onError:e=>{console.log("Error in clearMolsetWorkingCopy()",e)}})}}}),U={key:0,class:"fetching-err"},x={key:1,class:"fetching"},z=P({__name:"BaseFetching",props:{text:{default:"Fetching"},failText:{default:"Failed to fetch"},error:{type:Boolean,default:!1}},setup(e){const s=e;return(t,i)=>s.error?(d(),f("span",U,p(t.failText),1)):(d(),f("span",x,p(t.text),1))}}),T=b(z,[["__scopeId","data-v-9418e46f"]]);export{T as B,B as u};
