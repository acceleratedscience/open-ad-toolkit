import{d as I,a as D,c as S,r as _,w as V,f as w,A as b,k as v,i as d,l as h,u as x,h as C,j as y,x as t,v as a,B as p,t as F,n as B,p as M,q as E,F as R,C as P,D as z,E as A}from"./index-juGIK60X.js";import{_ as H}from"./JsonViewer.vue_vue_type_script_setup_true_lang-4fSEEgnR.js";import{M as K}from"./MolViewer-8oo8dSmH.js";import{B as N}from"./BaseFetching-muE6g1O5.js";import{i as T}from"./initRDKit-FDLbkFFb.js";import"./BreadCrumbs-5PNbBUT3.js";import"./TextBox-096FS_mU.js";import"./MolGridStore-YRzXPe7r.js";const $=I({__name:"TheMolFromIdentifier",props:{identifier:{}},setup(u){const c=D(),e=S(),n=u,s=_(!1),m=_("");n.identifier.startsWith("InChI=")?(e.setMolIdentifier("inchi",n.identifier),e.fetchMolVizData(n.identifier)):(i(n.identifier),e.fetchMolVizData(n.identifier)),f(n.identifier),V(()=>n.identifier,o=>{e.clear(),o&&f(o)});async function f(o=null){if(!o)return;let l=!1;return s.value=!0,m.value="",w(b.getMolData(o),{onSuccess:r=>{const g=!e.inchi;e.setMolData(r),l=!0,g&&e.fetchMolVizData(e.inchi)},loading:s,loadingError:m}),l}async function i(o){await T();let l=window.RDKit.get_mol(o);if(!l)return;const r=l.get_inchi();e.setMolIdentifier("canonical_smiles",o),e.setMolIdentifier("inchi",r),e.setMolIdentifier("inchikey",window.RDKit.get_inchikey_for_inchi(r))}return(o,l)=>v(e).loading?(d(),h(N,{key:0})):v(c).query.use=="json"?(d(),h(H,{key:1,data:v(e).mol},null,8,["data"])):(d(),h(K,{key:2,context:"identifier",loading:s.value,loadingError:m.value,onRetryLoad:l[0]||(l[0]=r=>f(o.identifier))},null,8,["loading","loadingError"]))}}),k=u=>(P("data-v-501187d9"),u=u(),z(),u),O=k(()=>t("h3",null,"Display any molecule",-1)),W=k(()=>t("br",null,null,-1)),Y={key:0,class:"error-msg"},j={class:"fields"},q=I({__name:"TheMolInput",setup(u){const c=x(),e=_(""),n=_("");function s(f){const i={inchi:"InChI=1S/C10H14O/c1-7-5-9(11)6-8(2)10(7,3)4/h5-6H,1-4H3",smiles:"CC1=CC(=O)C=C(C1(C)C)C",name:"penguinone",inchikey:"RHIYIMQPIGYWEK-UHFFFAOYSA-N",pid:"12564106"};e.value=i[f],n.value=""}function m(){e.value?c.push({name:"mol",params:{identifier:e.value.toString()}}):n.value="Please enter a molecule identifier."}return(f,i)=>{const o=C("cv-text-input"),l=C("cv-button");return d(),y(R,null,[O,t("p",null,[a(" Accepted identifiers are: "),t("b",null,[t("a",{href:"#",onClick:i[0]||(i[0]=p(r=>s("inchi"),["prevent"]))},"InChI")]),a(" or "),t("b",null,[t("a",{href:"#",onClick:i[1]||(i[1]=p(r=>s("smiles"),["prevent"]))},"SMILES")]),a("."),W,a(" When a molecule is listed on PubChem, you can also use its "),t("b",null,[t("a",{href:"#",onClick:i[2]||(i[2]=p(r=>s("name"),["prevent"]))},"name")]),a(", "),t("b",null,[t("a",{href:"#",onClick:i[3]||(i[3]=p(r=>s("inchikey"),["prevent"]))},"InChIKey")]),a(" or "),t("b",null,[t("a",{href:"#",onClick:i[4]||(i[4]=p(r=>s("pid"),["prevent"]))},"PID")]),a(". ")]),t("form",{id:"input-form",onSubmit:p(m,["prevent"])},[n.value?(d(),y("div",Y,F(n.value),1)):B("",!0),t("div",j,[M(o,{modelValue:e.value,"onUpdate:modelValue":i[5]||(i[5]=r=>e.value=r),type:"text",placeholder:"dopamine","hide-label":!0},null,8,["modelValue"]),M(l,{size:"default"},{default:E(()=>[a("Display")]),_:1})])],32)],64)}}}),L=A(q,[["__scopeId","data-v-501187d9"]]),te=I({__name:"MolPage",props:{identifier:{}},setup(u){return(c,e)=>c.identifier?(d(),h($,{key:0,identifier:c.identifier},null,8,["identifier"])):(d(),h(L,{key:1}))}});export{te as default};
