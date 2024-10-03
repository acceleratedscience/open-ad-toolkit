import{d as D,u as A,z as N,g as S,h as s,i as I,v as t,t as w,n as y,p as B,q as c,F as E,a as R,r as h,w as z,e as L,A as T,k,j as x,c as H,B as O,C as _,m as V,D as q,E as U,G as Y}from"./index-CoGQGb_u.js";import{u as j}from"./MolViewerStore-GYm8pbqh.js";import{_ as G}from"./JsonViewer.vue_vue_type_script_setup_true_lang-CEoPAOYf.js";import{M as Q}from"./MolViewer-BFbwHF1f.js";import{B as W}from"./BaseFetchingFile-BKvSjiJ2.js";import{i as X}from"./initRDKit-WaX80hiK.js";import"./BreadCrumbs-ULNMOtem.js";import"./TextBox-xPjm95_U.js";import"./BaseFetching-KBmQDor1.js";import"./MolRender3D--ExueadL.js";const J={class:"error-msg"},Z={class:"soft"},ee=t("br",null,null,-1),oe=D({__name:"BaseError",props:{loadingError:{},errorText:{default:"Something went wrong."},backLink:{}},setup(p){const d=A();N();const e=p;function l(){e.backLink?d.push(e.backLink):d.back()}return(i,v)=>{const m=S("cv-button");return s(),I(E,null,[t("div",J,w(i.errorText),1),t("p",Z,w(i.loadingError),1),ee,y(m,{onClick:l},{default:B(()=>[c("Back")]),_:1})],64)}}}),te=D({__name:"TheMolFromIdentifier",props:{identifier:{}},setup(p){const d=R(),e=j(),l=p,i=h(!1),v=h(!1),m=h(""),a=h(d.matched[0].name);l.identifier.startsWith("InChI=")?(e.setSmolIdentifier("inchi",l.identifier),e.fetchSmolVizData(l.identifier)):(o(l.identifier),e.fetchSmolVizData(l.identifier)),b(l.identifier),z(()=>l.identifier,n=>{e.clear(),n&&b(n)});async function b(n=null){n&&(i.value=!0,m.value="",a.value=a.value.replace(/^headless-/,""),a.value=="mol"?(v.value=!0,M(n,()=>{v.value=!1},()=>{g(n,()=>{v.value=!1},()=>{v.value=!1})})):["smol","mol"].includes(a.value)?M(n):a.value=="mmol"&&g(n))}function M(n,u=()=>{},f=()=>{}){L(T.getSmolData(n),{onSuccess:C=>{const r=!e.inchi;e.setMolData(C,"smol"),r&&e.fetchSmolVizData(e.inchi),u()},onError:f,loading:i,loadingError:m})}function g(n,u=()=>{},f=()=>{}){L(T.getMmolData(n),{onSuccess:C=>{e.setMolData(C,"mmol"),u()},onError:f,loading:i,loadingError:m})}async function o(n){await X();let u=window.RDKit.get_mol(n);if(!u)return;const f=u.get_inchi();e.setSmolIdentifier("canonical_smiles",n),e.setSmolIdentifier("inchi",f),e.setSmolIdentifier("inchikey",window.RDKit.get_inchikey_for_inchi(f))}return(n,u)=>m.value?(s(),k(oe,{key:0,loadingError:m.value,backLink:"/mol"},null,8,["loadingError"])):i.value||v.value?(s(),k(W,{key:1})):x(d).query.use=="json"?(s(),k(G,{key:2,data:x(e).smol},null,8,["data"])):(s(),k(Q,{key:3,context:"identifier",loading:i.value,loadingError:m.value,onRetryLoad:u[0]||(u[0]=f=>b(n.identifier))},null,8,["loading","loadingError"]))}}),F=p=>(q("data-v-8beda03f"),p=p(),U(),p),ne=F(()=>t("h3",null,"Display any molecule",-1)),le=F(()=>t("br",null,null,-1)),ie=F(()=>t("br",null,null,-1)),ae={key:0},re={key:1},se={key:0,class:"error-msg"},ue={class:"fields"},ce=D({__name:"TheMolInput",setup(p){const d=A(),e=h(""),l=h(""),i=h("smol"),v=[{label:"Small molecules",val:"smol"},{label:"Macromolecules",val:"mmol"}],m=H(()=>{switch(i.value){case"smol":return"InChI, SMILES, name, InChIKey or PubChem ID";case"mmol":return"FASTA sequence or PDB id";default:return"Molecule identifier"}});function a(g){const o={inchi:"InChI=1S/C10H14O/c1-7-5-9(11)6-8(2)10(7,3)4/h5-6H,1-4H3",smiles:"CC1=CC(=O)C=C(C1(C)C)C",name:"penguinone",inchikey:"RHIYIMQPIGYWEK-UHFFFAOYSA-N",cid:"681",fasta:"IINVKTSLKTIIKNALDKIQX",pdbId:"8t3n"};e.value=o[g],l.value=""}function b(){e.value="",l.value=""}function M(){e.value?d.push({name:i.value,params:{identifier:e.value.toString()}}):l.value="Please enter a molecule identifier."}return(g,o)=>{const n=S("cv-radio-button"),u=S("cv-radio-group"),f=S("cv-text-input"),C=S("cv-button");return s(),I(E,null,[ne,le,y(u,{legendText:"Search options",hideLegend:!0,onChange:b},{default:B(()=>[(s(),I(E,null,O(v,({label:r,val:K},P)=>y(n,{modelValue:i.value,"onUpdate:modelValue":o[0]||(o[0]=$=>i.value=$),key:P,name:"group-1",label:r,value:K,"hide-label":!1,"label-left":!1},null,8,["modelValue","label","value"])),64))]),_:1}),ie,["mol","smol"].includes(i.value)?(s(),I("p",ae,[c(" Accepted identifiers are: "),t("b",null,[t("a",{href:"#",onClick:o[1]||(o[1]=_(r=>a("inchi"),["prevent"]))},"InChI")]),c(" or "),t("b",null,[t("a",{href:"#",onClick:o[2]||(o[2]=_(r=>a("smiles"),["prevent"]))},"SMILES")]),c(" and "),t("b",null,[t("a",{href:"#",onClick:o[3]||(o[3]=_(r=>a("name"),["prevent"]))},"name")]),c(", "),t("b",null,[t("a",{href:"#",onClick:o[4]||(o[4]=_(r=>a("inchikey"),["prevent"]))},"InChIKey")]),c(" or "),t("b",null,[t("a",{href:"#",onClick:o[5]||(o[5]=_(r=>a("cid"),["prevent"]))},"PubChem CID")]),c(" when a molecule is listed on PubChem. ")])):V("",!0),["mol","mmol"].includes(i.value)?(s(),I("p",re,[c(" Macromolecules can be found by their "),t("b",null,[t("a",{href:"#",onClick:o[6]||(o[6]=_(r=>a("fasta"),["prevent"]))},"FASTA sequence")]),c(" or "),t("b",null,[t("a",{href:"#",onClick:o[7]||(o[7]=_(r=>a("pdbId"),["prevent"]))},"PDB id")]),c(". ")])):V("",!0),t("form",{id:"input-form",onSubmit:_(M,["prevent"])},[l.value?(s(),I("div",se,w(l.value),1)):V("",!0),t("div",ue,[y(f,{modelValue:e.value,"onUpdate:modelValue":o[8]||(o[8]=r=>e.value=r),type:"text",placeholder:m.value,"hide-label":!0},null,8,["modelValue","placeholder"]),y(C,{size:"default"},{default:B(()=>[c("Display")]),_:1})])],32)],64)}}}),de=Y(ce,[["__scopeId","data-v-8beda03f"]]),Ce=D({__name:"MolPage",props:{identifier:{}},setup(p){return(d,e)=>d.identifier?(s(),k(te,{key:0,identifier:d.identifier},null,8,["identifier"])):(s(),k(de,{key:1}))}});export{Ce as default};