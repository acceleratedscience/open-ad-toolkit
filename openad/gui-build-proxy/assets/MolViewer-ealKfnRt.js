import{d as E,g as ve,h as e,i as t,k as M,p as de,n as T,I as U,q as g,m as p,T as ft,S as O,C as ie,D as re,v as a,E as H,j as $,F as D,r as R,H as fe,P as dt,w as ut,Q as ne,Y as C,a1 as yt,A as F,t as y,u as _e,L as ye,b as J,a2 as gt,c as k,a3 as G,e as he,z as pe,a4 as kt,X as bt,R as vt,W as ge,B as wt,o as $t,a5 as St,a6 as Mt}from"./index-DwYRJoEC.js";import{u as j}from"./MolViewerStore-VkF0_V-U.js";import{u as ce,B as P}from"./BaseFetching-CtAuO-DI.js";import{B as Dt}from"./BreadCrumbs-CUWBa4am.js";import{M as _t}from"./MolRender3D-C-QWcCb6.js";const Ct=h=>(ie("data-v-bccbb3e2"),h=h(),re(),h),Bt=Ct(()=>a("div",{class:"filler"},null,-1)),It=E({__name:"BreadCrumbsNot",props:{backto:{},collapsed:{type:Boolean}},setup(h){return(d,i)=>{const s=ve("router-link");return e(),t("div",{id:"not-breadrumbs",class:O({collapsed:d.collapsed})},[d.backto?(e(),M(s,{key:0,to:d.backto,class:"backto"},{default:de(()=>[T(U,{icon:"icn-arrow-left"}),g(" Back")]),_:1},8,["to"])):p("",!0),Bt,ft(d.$slots,"default",{},void 0,!0)],2)}}}),Vt=H(It,[["__scopeId","data-v-bccbb3e2"]]),Tt=["innerHTML"],At=E({__name:"MolViewerVizSmol",setup(h){const d=j();return(i,s)=>(e(),t(D,null,[a("div",{class:"container-2d",innerHTML:$(d).svg},null,8,Tt),T(_t,{data3D:$(d).mdl,data3DFormat:"sdf",molName:$(d).name},null,8,["data3D","molName"])],64))}}),Ft=H(At,[["__scopeId","data-v-908ceb7a"]]),Et=E({__name:"MolViewerVizMmol",setup(h){const d=j();return(i,s)=>(e(),M(_t,{data3D:$(d).mmolData3D,data3DFormat:$(d).mmolData3DFormat,molName:$(d).name},null,8,["data3D","data3DFormat","molName"]))}}),Ht=H(Et,[["__scopeId","data-v-350c2e23"]]),Ot=h=>(ie("data-v-59dd251b"),h=h(),re(),h),Lt=Ot(()=>a("option",{value:"",disabled:"",hidden:""},null,-1)),Nt=["value","hidden"],xt=E({__name:"OverflowMenu",props:{options:{},disabled:{type:Boolean}},setup(h){const d=h,i=R(""),s=R(!1);return fe(async()=>{const l=d.options.find(o=>o.selected);l&&(i.value=l.val),await dt(),s.value=!0}),ut(()=>i.value,async l=>{if(s.value){if(l){const o=d.options.find(u=>u.val===l);if(o){const u=o.action;u()}}await dt(),i.value=""}}),(l,o)=>(e(),t("div",{class:O(["btn-wrap",{disabled:l.disabled}])},[T(ne,{icon:"icn-overflow"}),C(a("select",{"onUpdate:modelValue":o[0]||(o[0]=u=>i.value=u)},[Lt,(e(!0),t(D,null,F(l.options,(u,m)=>(e(),t("option",{key:m,value:u.val,hidden:!!u.hide},y(u.disp),9,Nt))),128))],512),[[yt,i.value]])],2))}}),zt=H(xt,[["__scopeId","data-v-59dd251b"]]);function mt(){const h=_e(),d=ye(),i=J(),s=ce(),l=j();return function(m,v=!1,c={}){const n=(i.moduleName=="MolViewer"?i.filenameNaked:c.defaultName)||"untitled";return d.display("ModalSaveFile",{path:i.pathDir,filename:n,dataType:m,exportOptions:v},{onSubmit:o})};async function o({destinationPath:u,ext:m,srcDataType:v},c=!1){try{v=="molset"?m=="molset.json"?await s.saveMolsetAsJSON(u,{newFile:!0,force:c}):m=="sdf"?await s.saveMolsetAsSDF(u,{newFile:!0,force:c}):m=="csv"?await s.saveMolsetAsCSV(u,{newFile:!0,force:c}):m=="smi"&&await s.saveMolsetAsSmiles(u,{newFile:!0,force:c}):v=="smol"?m=="smol.json"?await l.saveSmolAsJSON(u,{newFile:!0,force:c}):m=="sdf"?await l.saveSmolAsSDF(u,{newFile:!0,force:c}):m=="csv"?await l.saveSmolAsCSV(u,{newFile:!0,force:c}):m=="mol"?await l.saveSmolAsMDL(u,{newFile:!0,force:c}):m=="smi"&&await l.saveSmolAsSMILES(u,{newFile:!0,force:c}):["cif","pdb"].includes(v)&&(m=="mmol.json"?await l.saveMmolAsMmolJson(u,{newFile:!0,force:c}):m=="cif"?await l.saveMmolAsCIF(u,{newFile:!0,force:c}):m=="pdb"&&await l.saveMmolAsPDB(u,{newFile:!0,force:c}));const n=gt(u);h.push("/~/"+n)}catch(n){const r=(n==null?void 0:n.data)||"An error occurred while saving the file.",b=n==null?void 0:n.status;n!=null&&n.status||console.error(n),b==409?d.confirm(r,{title:"Overwrite existing file?",onSubmit:()=>o({destinationPath:u,ext:m,srcDataType:v},!0)}):d.alert(r,{title:"Error"})}}}const Rt=E({__name:"OverflowMenuMol",props:{dataType:{},disabled:{type:Boolean}},setup(h){const d=_e(),i=J(),s=j(),l=ye(),o=mt(),u=h,m=[{val:"save-as",disp:"Save as...",action:v}];i.active&&["smol","mmol","cif","pdb"].includes(i.fileType)&&m.push({val:"delete",disp:"Delete",action:c});function v(){if(console.log("dataType: ",u.dataType),u.dataType=="smol")return o("smol",!0,{defaultName:s.nameSlug});if(u.dataType=="cif")return o("cif",!0,{defaultName:s.nameSlug});if(u.dataType=="pdb")return o("pdb",!0,{defaultName:s.nameSlug});console.error(`actionSaveAs() -> props.dataType "${u.dataType}" not recognized`)}async function c(){l.confirm("This cannot be undone.",{title:"Please confirm deletion",primaryBtn:"Delete",onSubmit:async()=>{await i.deleteFile(i.pathAbsolute)&&d.push("/~/"+i.pathDir)}})}return(n,r)=>(e(),M(zt,{options:m,disabled:n.disabled},null,8,["disabled"]))}}),jt=H(Rt,[["__scopeId","data-v-67ccb099"]]),Pt=E({__name:"BaseIconButtonToggle",props:{icon:{},iconHover:{},iconOn:{},iconOnHover:{},color:{},colorHover:{},colorOn:{},colorOnHover:{},mini:{type:Boolean}},emits:["toggle-on","toggle-off"],setup(h,{expose:d,emit:i}){const s=i,l=h,o=R(!1),u=R(!1),m={soft:"rgba(0,0,0,.3)",semiSoft:"rgba(0,0,0,.6)",hard:"#393939"},v=k(()=>{const n={};return n["--btn-color"]=l.color?l.color:m.soft,n["--btn-color-hover"]=l.colorHover?l.colorHover:m.soft,n["--btn-color-on"]=l.colorOn?l.colorOn:m.hard,n["--btn-color-on-hover"]=l.colorOnHover?l.colorOnHover:m.hard,n});function c(n=void 0,r=!1){n!==void 0?o.value=n:o.value=!o.value,o.value?(u.value=!0,r||s("toggle-on")):r||s("toggle-off")}return d({toggle:c}),(n,r)=>(e(),t("div",{class:O(["icn-btn",{[l.icon]:!0,mini:l.mini,"toggle-on":o.value,"ignore-hover":u.value,"has-hover-icn":l.iconHover,"has-on-icn":l.iconOn,"has-on-hover-icn":l.iconOnHover,"has-custom-color":!!l.color,"has-custom-hover-color":!!l.colorHover,"has-custom-toggle-color":!!l.colorOn}]),onClick:r[0]||(r[0]=b=>c()),onMouseleave:r[1]||(r[1]=b=>u.value=!1),style:G(v.value)},[T(U,{class:"base-icn",icon:l.icon},null,8,["icon"]),l.iconHover?(e(),M(U,{key:0,class:"hover-icn",icon:l.iconHover},null,8,["icon"])):p("",!0),l.iconOn?(e(),M(U,{key:1,class:"on-icn",icon:l.iconOn},null,8,["icon"])):p("",!0),l.iconOnHover?(e(),M(U,{key:2,class:"on-hover-icn",icon:l.iconOnHover},null,8,["icon"])):p("",!0)],38))}}),Wt=H(Pt,[["__scopeId","data-v-2dfe3d6f"]]),Ut=E({__name:"BaseBookmark",props:{mol:{}},setup(h){const d=h,i=R(null),s=R(!1),l=R("");ut(()=>d.mol,m),fe(m);async function o(){if(!d.mol){l.value="x";return}he(pe.addMolToMyMols(d.mol),{onError:()=>{i.value&&i.value.toggle(!1,!0)},loading:s,loadingError:l})}async function u(){if(!d.mol){l.value="x";return}he(pe.removeMolFromMyMols(d.mol),{onError:v=>{i.value&&i.value.toggle(v.status,!0)},loading:s,loadingError:l})}function m(){!d.mol||!("properties"in d.mol)||d.mol&&he(pe.checkMolInMyMols(d.mol),{onSuccess:v=>{i.value&&i.value.toggle(v.status,!0)}})}return(v,c)=>(e(),M(Wt,{ref_key:"toggleBtn",ref:i,class:O({loading:s.value,error:!!l.value}),icon:"icn-bookmark",iconHover:"icn-bookmark-full",iconOn:"icn-bookmark-full",iconOnHover:"icn-bookmark-remove",onToggleOn:o,onToggleOff:u},null,8,["class"]))}}),Kt=H(Ut,[["__scopeId","data-v-c39285cc"]]),Gt={class:"display"},Jt=["value"],qt={key:0,value:"before",disabled:""},Qt=["value"],Xt={key:1,value:"after",disabled:""},Yt=E({__name:"BasePagination",props:{total:{},modelValue:{default:1},disabled:{type:Boolean}},emits:["update:modelValue"],setup(h,{emit:d}){const i=d,s=h,l=k(()=>kt(s.total)),o=k(()=>{const n=[],r=Math.max(1,s.modelValue-100),b=Math.min(s.total,r+200);for(let w=r;w<=b;w++)n.push(w);return n});function u(){s.modelValue>1&&i("update:modelValue",s.modelValue-1)}function m(){s.modelValue<s.total&&i("update:modelValue",s.modelValue+1)}function v(c){c&&(c=typeof c=="string"?parseInt(c):c,i("update:modelValue",c))}return(c,n)=>(e(),t("div",{class:O(["pagination",{disabled:c.disabled}])},[a("div",{class:"btn",onClick:u},[T(U,{icon:"icn-caret-left"})]),a("div",{class:"btn",onClick:m},[T(U,{icon:"icn-caret-right"})]),a("div",Gt,[a("span",null,y(s.modelValue)+" / "+y(l.value),1),a("select",{value:s.modelValue,onChange:n[0]||(n[0]=r=>{var b;return v((b=r.target)==null?void 0:b.value)})},[o.value[0]>1?(e(),t("option",qt,"...")):p("",!0),(e(!0),t(D,null,F(o.value,r=>(e(),t("option",{key:r,value:r},y(r),9,Qt))),128)),o.value[o.value.length-1]<s.total?(e(),t("option",Xt,"...")):p("",!0)],40,Jt)])],2))}}),Zt=H(Yt,[["__scopeId","data-v-fa8eea59"]]),eo=h=>(ie("data-v-d0ded2e3"),h=h(),re(),h),to={id:"title-wrap"},oo={class:"v-align"},so=eo(()=>a("div",{class:"filler"},null,-1)),ao=E({__name:"MolViewerTitle",props:{context:{default:"file"},loading:{type:Boolean}},setup(h){const d=_e(),i=J(),s=j(),l=ce(),o=h,u=k(()=>o.loading?"Loading":s.name),m=k(()=>o.context=="molset"?"icn-file-molset":s.isSmol?"icn-file-smol":"icn-file-mmol"),v=k({get:()=>s.molFromMolsetIndex||1,set:s.setMolFromMolsetIndex});return(c,n)=>(e(),t("div",to,[a("div",oo,[T(U,{class:O(["icn-mol",{loading:c.loading}]),icon:m.value,size:"large"},null,8,["class","icon"])]),a("h2",{id:"data-name","data-val":"{{ molName }}",class:O({loading:c.loading})},y($(bt)(u.value)),3),so,$(s).smol?(e(),M(Kt,{key:0,mol:$(s).smol},null,8,["mol"])):p("",!0),T(jt,{dataType:$(s).dataType,disabled:!!c.loading},null,8,["dataType","disabled"]),c.context=="molset"?(e(),M(Zt,{key:1,modelValue:v.value,"onUpdate:modelValue":n[0]||(n[0]=r=>v.value=r),total:$(l).total},null,8,["modelValue","total"])):p("",!0),c.context=="molset"?(e(),M(ne,{key:2,icon:"icn-close",icnSize:"small",btnStyle:"carbon",onClick:n[1]||(n[1]=r=>$(s).setMolFromMolsetIndex(null))})):(e(),M(ne,{key:3,icon:"icn-close",icnSize:"small",btnStyle:"default",onClick:n[2]||(n[2]=r=>c.context=="file"?$(i).exitViewer():c.context=="identifier"?$(d).push({name:"mol"}):null),style:{"margin-right":"-8px"}}))]))}}),lo=H(ao,[["__scopeId","data-v-d0ded2e3"]]),no=E({__name:"TheButtonSaveMol",setup(h){const d=j(),i=J(),s=ce(),l=mt(),o=R(!1),u=k(()=>i.active&&!["smol","molset"].includes(i.fileType||"")),m=k(()=>{if(s.active){if(s.context=="json"||s.context=="my-mols")return!1}else if(i.fileType=="smol")return!1;return!0}),v=k(()=>d.hasChanges);async function c(){let n=!1;o.value=!0,s.active?s.context=="json"?n=await s.replaceMolInMolset(i.path,d.smol,s.context):s.context=="my-mols"?n=await s.replaceMolInMolset(i.path,d.smol,s.context):(s.setHasChanges(!1),n=await l("smol",!0,{defaultName:d.nameSlug}),s.setHasChanges(!n)):!i.active||u.value?(d.setHasChanges(!1),n=await l("smol",!0,{defaultName:d.nameSlug}),d.setHasChanges(!n)):i.fileType=="smol"&&(n=await d.saveSmolAsJSON(i.path,{newFile:!1})),n&&(d.setHasChanges(!1),s.setHasChanges(!1)),o.value=!1}return(n,r)=>{const b=ve("cv-button");return v.value?(e(),M(b,{key:0,size:"small",kind:"primary",class:O(["btn-save",{saving:o.value,"save-as":m.value}]),disabled:o.value,onClick:c},null,8,["class","disabled"])):p("",!0)}}}),io=H(no,[["__scopeId","data-v-07df3337"]]),ro=E({__name:"TheButtonEnrichMol",setup(h){const d=j(),i=ce(),s=R(!1),l=k(()=>!d.enriched);async function o(){s.value=!0,await d.enrichSmol()&&d.molFromMolset&&i.setHasChanges(!0),s.value=!1}return(u,m)=>{const v=ve("cv-button");return l.value?(e(),M(v,{key:0,size:"small",kind:"secondary",class:O(["btn-enrich",{enriching:s.value}]),disabled:s.value,title:"Calculate molecular properties and load data from PubChem",onClick:o},null,8,["class","disabled"])):p("",!0)}}}),co=H(ro,[["__scopeId","data-v-7283d0a8"]]),K=h=>(ie("data-v-1fe8a102"),h=h(),re(),h),uo={id:"identification"},vo={key:0},_o=["data-copy"],mo={key:0,id:"data-inchi"},ho={key:1,class:"blank"},po={key:1},fo=["data-copy"],yo={key:0,id:"data-inchikey"},go={key:1,class:"blank"},ko={key:2},bo=["data-copy"],wo={id:"data-smiles"},$o={key:3},So=["data-copy"],Mo={key:0,id:"data-canonical-smiles"},Do={key:1,class:"blank"},Co={key:4},Bo=["data-copy"],Io={key:0,id:"data-isomeric-smiles"},Vo={key:1,class:"blank"},To={key:5},Ao=["data-copy"],Fo={key:0,id:"data-isomeric-smiles"},Eo={key:1,class:"blank"},Ho={key:6},Oo=["data-copy"],Lo=["href"],No={key:1,class:"blank"},xo=K(()=>a("hr",null,null,-1)),zo={key:1,id:"fetch-fail"},Ro=K(()=>a("div",{class:"error-msg"},"Something went wrong fetching the molecule data.",-1)),jo={key:0,class:"status-msg"},Po=K(()=>a("h3",null,"Synonyms",-1)),Wo={class:"flip-v"},Uo={key:1,class:"cloak"},Ko=["title"],Go=K(()=>a("hr",null,null,-1)),Jo={id:"properties"},qo=K(()=>a("h3",null,"Properties",-1)),Qo=["title"],Xo=["data-copy"],Yo=K(()=>a("div",{class:"filler"},null,-1)),Zo={class:"val"},es={key:0},ts={key:1,id:"analysis"},os=K(()=>a("h3",null,"Analysis",-1)),ss=["data-copy"],as=K(()=>a("br",null,null,-1)),ls=250,ns=150,is=E({__name:"MolViewerDataSmol",props:{loading:{type:Boolean},loadingError:{}},emits:["retryLoad"],setup(h,{emit:d}){const i=j(),s=vt(),l=d,o=k(()=>i.smol),u=k(()=>!!(i.smol&&"analysis"in i.smol&&i.smol.analysis.length)),m=k(()=>{var _;return"synonyms"in o.value&&(_=o.value)!=null&&_.synonyms?o.value.synonyms.length:0}),v=k(()=>`${Math.ceil(m.value/c.value)*22}px`),c=k(()=>s.contentWidth?Math.floor(s.contentWidth/ns)||1:4),n=k(()=>`calc((100% - ${(c.value-1)*20}px) / ${c.value})`),r=k(()=>{const _=5*c.value+c.value;return m.value>_}),b=k(()=>r.value?"110px":""),w=k(()=>s.contentWidth?Math.floor(s.contentWidth/ls)||1:3),z=k(()=>`calc((100% - ${(w.value-1)*40}px) / ${w.value})`),f=k(()=>{var I;if(!("properties"in o.value))return{};const _=(I=o.value)!=null&&I.properties?Object.keys(o.value.properties).length:0,A=Math.ceil(_/w.value)*22;return A?{height:`${A}px`}:{}});function S(_){_.currentTarget.classList.toggle("expand")}return(_,A)=>{var q,Q,X,Y,Z,ee,te,oe,se,ae,le,V,N,ke,be,we,$e,Se,Me,De,Ce,Be,Ie,Ve,Te,Ae,Fe,Ee,He,Oe,Le,Ne,xe,ze,Re,je,Pe,We,Ue,Ke,Ge,Je,qe,Qe,Xe,Ye,Ze,et,tt,ot,st,at,lt,nt,it,rt;const I=ve("cv-button"),B=ge("click-to-copy");return e(),t(D,null,[a("div",uo,[$(i).enriched||(Q=(q=o.value)==null?void 0:q.identifiers)!=null&&Q.inchi?(e(),t("div",vo,[C((e(),t("b",{"data-copy":`InChI: ${(Y=(X=o.value)==null?void 0:X.identifiers)==null?void 0:Y.inchi}`},[g("InChI: ")],8,_o)),[[B]]),(ee=(Z=o.value)==null?void 0:Z.identifiers)!=null&&ee.inchi?C((e(),t("span",mo,[g(y((oe=(te=o.value)==null?void 0:te.identifiers)==null?void 0:oe.inchi),1)])),[[B]]):(e(),t("span",ho,"-")),_.loading?(e(),M(P,{key:2,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),$(i).enriched||(ae=(se=o.value)==null?void 0:se.identifiers)!=null&&ae.inchikey?(e(),t("div",po,[C((e(),t("b",{"data-copy":`InChIKey: ${(V=(le=o.value)==null?void 0:le.identifiers)==null?void 0:V.inchikey}`},[g("InChIKey: ")],8,fo)),[[B]]),(ke=(N=o.value)==null?void 0:N.identifiers)!=null&&ke.inchikey?C((e(),t("span",yo,[g(y((we=(be=o.value)==null?void 0:be.identifiers)==null?void 0:we.inchikey),1)])),[[B]]):(e(),t("span",go,"-")),_.loading?(e(),M(P,{key:2,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),(Se=($e=o.value)==null?void 0:$e.identifiers)!=null&&Se.smiles&&(!((De=(Me=o.value)==null?void 0:Me.identifiers)!=null&&De.canonical_smiles)||(Be=(Ce=o.value)==null?void 0:Ce.identifiers)!=null&&Be.isomeric_smiles)?(e(),t("div",ko,[C((e(),t("b",{"data-copy":`SMILES: ${(Ve=(Ie=o.value)==null?void 0:Ie.identifiers)==null?void 0:Ve.smiles}`},[g("SMILES: ")],8,bo)),[[B]]),C((e(),t("span",wo,[g(y((Ae=(Te=o.value)==null?void 0:Te.identifiers)==null?void 0:Ae.smiles),1)])),[[B]]),_.loading?(e(),M(P,{key:0,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),$(i).enriched||(Ee=(Fe=o.value)==null?void 0:Fe.identifiers)!=null&&Ee.canonical_smiles?(e(),t("div",$o,[C((e(),t("b",{"data-copy":`Canonical SMILES: ${(Oe=(He=o.value)==null?void 0:He.identifiers)==null?void 0:Oe.canonical_smiles}`},[g("Canonical SMILES: ")],8,So)),[[B]]),(Ne=(Le=o.value)==null?void 0:Le.identifiers)!=null&&Ne.canonical_smiles?C((e(),t("span",Mo,[g(y((ze=(xe=o.value)==null?void 0:xe.identifiers)==null?void 0:ze.canonical_smiles),1)])),[[B]]):(e(),t("span",Do,"-")),_.loading?(e(),M(P,{key:2,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),$(i).enriched||(je=(Re=o.value)==null?void 0:Re.identifiers)!=null&&je.isomeric_smiles?(e(),t("div",Co,[C((e(),t("b",{"data-copy":`Isomeric SMILES: ${(We=(Pe=o.value)==null?void 0:Pe.identifiers)==null?void 0:We.isomeric_smiles}`},[g("Isomeric SMILES: ")],8,Bo)),[[B]]),(Ke=(Ue=o.value)==null?void 0:Ue.identifiers)!=null&&Ke.isomeric_smiles?C((e(),t("span",Io,[g(y((Je=(Ge=o.value)==null?void 0:Ge.identifiers)==null?void 0:Je.isomeric_smiles),1)])),[[B]]):(e(),t("span",Vo,"-")),_.loading?(e(),M(P,{key:2,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),$(i).enriched||(Qe=(qe=o.value)==null?void 0:qe.identifiers)!=null&&Qe.molecular_formula?(e(),t("div",To,[C((e(),t("b",{"data-copy":`Formula: ${(Ye=(Xe=o.value)==null?void 0:Xe.identifiers)==null?void 0:Ye.molecular_formula}`},[g("Formula: ")],8,Ao)),[[B]]),(et=(Ze=o.value)==null?void 0:Ze.identifiers)!=null&&et.molecular_formula?C((e(),t("span",Fo,[g(y((ot=(tt=o.value)==null?void 0:tt.identifiers)==null?void 0:ot.molecular_formula),1)])),[[B]]):(e(),t("span",Eo,"-")),_.loading?(e(),M(P,{key:2,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),$(i).enriched||(st=o.value.identifiers)!=null&&st.cid?(e(),t("div",Ho,[C((e(),t("b",{"data-copy":`PubChem CID: ${o.value.identifiers.cid}`},[g("PubChem CID: ")],8,Oo)),[[B]]),(lt=(at=o.value)==null?void 0:at.identifiers)!=null&&lt.cid?(e(),t("a",{key:0,id:"data-cid",href:`https://pubchem.ncbi.nlm.nih.gov/compound/${(nt=o.value.identifiers)==null?void 0:nt.cid}`,target:"_blank"},y(o.value.identifiers.cid),9,Lo)):(e(),t("span",No,"-")),_.loading?(e(),M(P,{key:2,text:"",failText:"x",error:!!_.loadingError},null,8,["error"])):p("",!0)])):p("",!0),!_.loading&&!_.loadingError?(e(),t(D,{key:7},[$(i).isSmol?(e(),M(co,{key:0})):p("",!0),T(io)],64)):p("",!0)]),xo,_.loading&&!_.loadingError?(e(),M(P,{key:0,text:"Fetching molecule data",error:!!_.loadingError},null,8,["error"])):_.loadingError?(e(),t("div",zo,[Ro,_.loadingError?(e(),t("div",jo,y(_.loadingError),1)):p("",!0),a("div",null,[T(I,{kind:"danger",size:"field",onClick:A[0]||(A[0]=L=>l("retryLoad"))},{default:de(()=>[g("Retry")]),_:1})])])):(e(),t(D,{key:2},[a("div",{id:"synonyms",style:G({"--truncated-height":b.value})},[Po,a("div",Wo,[r.value?(e(),t("a",{key:0,href:"#",class:"toggle-expand",onClick:wt(S,["prevent"])})):p("",!0),m.value&&"synonyms"in o.value?(e(),t("div",Uo,[a("div",{class:"synonyms-wrap",style:G({height:v.value})},[(e(!0),t(D,null,F((it=o.value)==null?void 0:it.synonyms,(L,W)=>C((e(),t("div",{key:W,title:L,style:G({width:n.value})},[g(y(L),1)],12,Ko)),[[B]])),128))],4)])):_.loading?(e(),M(P,{key:2,error:!!_.loadingError},null,8,["error"])):(e(),t(D,{key:3},[g("No synonyms available")],64))])],4),Go,a("div",Jo,[qo,"properties"in o.value?(e(),t("div",{key:0,class:"param-wrap",style:G(f.value)},[(e(!0),t(D,null,F((rt=o.value)==null?void 0:rt.properties,(L,W)=>(e(),t("div",{key:W,title:$(i).propertiesString?$(i).propertiesString[W]:"",class:O({empty:!L&&L!==0}),style:G({width:z.value})},[C((e(),t("div",{"data-copy":`${W}: ${L}`,class:"key"},[g(y(W)+":",1)],8,Xo)),[[B]]),Yo,C((e(),t("div",Zo,[g(y(L||L===0?L:"-"),1)])),[[B]])],14,Qo))),128))],4)):p("",!0)]),u.value?(e(),t("hr",es)):p("",!0),u.value&&$(i).smol&&"analysis"in $(i).smol?(e(),t("div",ts,[os,(e(!0),t(D,null,F($(i).smol.analysis,(L,W)=>(e(),t("div",{key:W,class:"item"},[a("details",null,[a("summary",null,[a("b",null,y(L.toolkit),1),g(" / "+y(L.function),1)]),(e(!0),t(D,null,F(L.results,(ht,pt)=>(e(),t("div",{key:pt,class:"result"},[(e(!0),t(D,null,F(ht,(ct,me)=>(e(),t("div",{key:me},[C((e(),t("b",{"data-copy":`${me}: ${ct}`},[g(y(me)+":",1)],8,ss)),[[B]]),g(),C((e(),t("span",null,[g(y(ct),1)])),[[B]]),as]))),128))]))),128))])]))),128))])):p("",!0)],64))],64)}}}),rs=H(is,[["__scopeId","data-v-1fe8a102"]]),cs=["data-copy"],ds={key:0},us=["title"],vs=["href"],_s=["data-copy"],ms=E({__name:"TableData",props:{data:{},inline:{type:Boolean,default:!1},header:{type:Boolean,default:!0},allowCopy:{type:Boolean,default:!1}},setup(h){const d=R(null),i=h,s=k(()=>{if(typeof i.data=="object"&&i.data!==null){let v=0,c=!0;for(const n in i.data){if(+n!==v){c=!1;break}v++}return c?"A":"B"}return null}),l=k(()=>m(i.data,s.value)),o=k(()=>{let v="",c=0;for(const n in l.value)(c>0||i.header)&&(v+=l.value[n].join(",")+`
`),c++;return v}),u=k(()=>{const v={header:[],body:[]};for(let c=0;c<l.value.length;c++){const n=l.value[c];Array.isArray(n)&&(c===0?n.forEach(r=>{v.header.push(r)}):v.body.push(n))}return v});function m(v,c){if(!v)return console.error("restructureData() - No data provided"),[];if(!c)return console.error("restructureData() - No data structure provided"),[];const n=[];if(c=="A"){let r=0;for(const b in v){const w=v[b];r===0&&(n[0]=[]),n[r+1]=[];for(const z in w)r===0&&n[r].push(z),n[r+1].push(w[z]);r++}}else if(c=="B"){let r=0;for(const b in v)r===0&&(n[0]=[]),n[r+1]=[b,v[b]],r++}return n}return(v,c)=>{const n=ge("click-to-copy");return e(),t("div",{ref_key:"overflowWrap",ref:d,class:O(["table-overflow-wrap",{inline:v.inline}])},[C((e(),t("table",{"data-copy":v.allowCopy?o.value:null,class:O({"key-val":s.value=="B"})},[v.header?(e(),t("thead",ds,[a("tr",null,[(e(!0),t(D,null,F(u.value.header,(r,b)=>(e(),t("th",{key:b},[a("span",{class:O({soft:!r})},y(r||"-"),3)]))),128))])])):p("",!0),a("tbody",null,[(e(!0),t(D,null,F(u.value.body,(r,b)=>(e(),t("tr",{key:b},[(e(!0),t(D,null,F(r,(w,z)=>(e(),t("td",{key:z,title:s.value=="B"&&String(w).length>20?String(w):""},[w&&typeof w=="string"&&w.match(/^http(s)?:\/\//)?(e(),t("a",{key:0,href:w,target:"_blank"},y(w),9,vs)):C((e(),t("span",{key:1,"data-copy":z==0?`${w}: ${r[1]}`:z==1?w:null,class:O({soft:!w})},[g(y(w||"-"),1)],10,_s)),[[n,s.value=="B"]])],8,us))),128))]))),128))])],10,cs)),[[n,v.allowCopy]])],2)}}}),ue=H(ms,[["__scopeId","data-v-0a093fa5"]]),hs={key:0},ps={key:1},fs={class:"small soft"},ys={class:"small soft vector"},gs=E({__name:"MolViewerDataMmolMatrices",props:{data:{}},setup(h){const d=h,i=k(()=>{const s=d.data;return Array.isArray(s)?s:[s]});return(s,l)=>(e(!0),t(D,null,F(i.value,(o,u)=>(e(),t("div",{key:u,class:"matrix-wrap-wrap"},[(e(!0),t(D,null,F(o.matrices,(m,v)=>{var c;return e(),t("div",{key:v,class:"matrix-wrap"},[(c=o.fields)!=null&&c.id?(e(),t("h5",hs,"Matrix "+y(o.fields.id),1)):p("",!0),(e(!0),t(D,null,F(o.fields,(n,r)=>(e(),t("div",{key:r},[g(" • "),a("b",null,y(r)+":",1),g(" "+y(n),1)]))),128)),Object.keys(o.fields).length?(e(),t("br",ps)):p("",!0),a("div",fs,[g(" Matrix"),String(v)!="_"?(e(),t(D,{key:0},[g(" - "+y(v),1)],64)):p("",!0)]),T(ue,{data:m,allowCopy:!0,header:!1},null,8,["data"]),a("div",ys,[g(" Vector"),String(v)!="_"?(e(),t(D,{key:0},[g(" - "+y(v),1)],64)):p("",!0)]),T(ue,{data:[o.vectors[v]],allowCopy:!0,header:!1},null,8,["data"])])}),128))]))),128))}}),ks=H(gs,[["__scopeId","data-v-aa496338"]]),x=h=>(ie("data-v-5c5b6fff"),h=h(),re(),h),bs=["innerHTML"],ws=x(()=>a("hr",null,null,-1)),$s={class:"data-block"},Ss=x(()=>a("h3",null,"Overview",-1)),Ms={class:"key-val"},Ds=["data-copy"],Cs=["href"],Bs={key:1},Is={class:"key-val"},Vs=["data-copy"],Ts=["href"],As={key:1},Fs={class:"key-val"},Es=x(()=>a("div",null,"FASTA",-1)),Hs=["href"],Os=["href"],Ls={key:1},Ns={class:"key-val"},xs=["data-copy"],zs={class:"key-val"},Rs=["data-copy"],js={class:"key-val inline"},Ps=["data-copy"],Ws={key:0},Us=["href"],Ks={key:0},Gs={key:1},Js={class:"key-val inline"},qs=["data-copy"],Qs={key:0},Xs=["href"],Ys={key:0},Zs={key:1},ea={key:0,class:"key-val"},ta=["data-copy"],oa=x(()=>a("hr",null,null,-1)),sa={class:"data-block"},aa=x(()=>a("h3",null,"Experimental Data Snapshot",-1)),la={class:"key-val"},na=x(()=>a("div",null,"Method",-1)),ia={class:"key-val"},ra=x(()=>a("div",null,"Resolution",-1)),ca={class:"key-val"},da=x(()=>a("div",null,"R-Value Free",-1)),ua={class:"key-val"},va=x(()=>a("div",null,"R-Value Work",-1)),_a={class:"key-val"},ma=x(()=>a("div",null,"R-Value Observed",-1)),ha=x(()=>a("hr",null,null,-1)),pa={key:0,class:"data-block"},fa=x(()=>a("h3",null,"Validation Report",-1)),ya=["src"],ga={key:1},ka={class:"data-block"},ba=x(()=>a("h3",null,"Meta Data",-1)),wa={key:0},$a=["href","name"],Sa={key:1},Ma={key:2},Da=E({__name:"MolViewerDataMmol",setup(h){const d=j(),i=J(),s=k(()=>i.ext=="json"&&i.ext2=="mmol"),l=k(()=>d.mmolData),o=k(()=>d.mmolDataHuman),u=k(()=>{var f,S;return((S=(f=l.value)==null?void 0:f.entry)==null?void 0:S.id)||null}),m=k(()=>{var S;let f=(S=l.value)!=null&&S.database_2?o.value["Database 2"]:null;if(f){for(let _=0;_<f.length;_++)if(f[_]["Database ID"]=="Pdb")return f[_].Doi}return null}),v=k(()=>{var S,_;let f=((_=(S=l.value)==null?void 0:S.pdbx_database_status)==null?void 0:_.recvd_initial_deposition_date)||null;return f?new Date(f).toLocaleDateString("en-US",{month:"short",day:"numeric",year:"numeric"}):null}),c=k(()=>{var A;const f=((A=l.value)==null?void 0:A.pdbx_audit_revision_history)||null;if(!f)return null;const S=Array.isArray(f)?f:[f];let _=S?S[0].revision_date:null;return _?new Date(_).toLocaleDateString("en-US",{month:"short",day:"numeric",year:"numeric"}):null}),n=k(()=>{var _;const f=(_=l.value)!=null&&_.audit_author?o.value["Audit Author"]:null;if(!f)return[];let S=Array.isArray(f)?f:[f];return S=S?S.map(A=>{let I=A.Name;return I=I.split(",").reverse(),I=I.map(B=>B.trim()),I[0].match(/^([a-zA-Z]\.)+$/)&&(I[0]=I[0].toUpperCase()),I=I.join(" "),I}):[],S}),r=k(()=>{var S,_;let f=(_=(S=l.value)==null?void 0:S.struct_keywords)!=null&&_.text?o.value["Structure Keywords"].Text:null;return f&&(f=f.split(","),f=f.map(A=>A.trim())),f}),b=k(()=>{if(!u.value)return null;const f=u.value.toLowerCase();return`https://files.rcsb.org/pub/pdb/validation_reports/${f.slice(1,3)}/${f}/${f}_multipercentile_validation.png`});fe(()=>{const f=window.location.hash;window.location.hash="",setTimeout(()=>{w(f)},1)});function w(f){if(!f)return;window.location.hash=f;const S=document.querySelector(`a[href="${f}"]`),_=S==null?void 0:S.closest("details");_&&_.setAttribute("open","true"),S==null||S.blur()}function z(f){return f?"https://scholar.google.com/scholar?q="+f.replace(/\s/g,"+"):""}return(f,S)=>{var A,I,B,q,Q,X,Y,Z,ee,te,oe,se,ae,le;const _=ge("click-to-copy");return e(),t(D,null,[a("div",{class:"capitalize",innerHTML:((I=(A=o.value)==null?void 0:A.Structure)==null?void 0:I.Title)||"<div class='soft'>No description available</div>"},null,8,bs),ws,a("div",$s,[Ss,a("div",Ms,[C((e(),t("div",{"data-copy":`PDB Entry: https://www.rcsb.org/structure/${u.value}`},[g("PDB Entry")],8,Ds)),[[_]]),u.value?(e(),t("a",{key:0,href:`https://www.rcsb.org/structure/${u.value}`,target:"_blank"},y(u.value),9,Cs)):(e(),t("div",Bs,"-"))]),a("div",Is,[C((e(),t("div",{"data-copy":`DOI Entry: https://doi.org/${m.value}`},[g("DOI")],8,Vs)),[[_]]),m.value?(e(),t("a",{key:0,href:`https://doi.org/${m.value}`,target:"_blank"},y(m.value),9,Ts)):(e(),t("div",As,"-"))]),a("div",Fs,[Es,u.value?(e(),t(D,{key:0},[a("a",{href:`https://www.rcsb.org/fasta/entry/${u.value}/display`,target:"_blank"},"View",8,Hs),g(" /  "),a("a",{href:`https://www.rcsb.org/fasta/entry/${u.value}`},"Download",8,Os)],64)):(e(),t("div",Ls,"-"))]),a("div",Ns,[C((e(),t("div",{"data-copy":`Deposited: ${v.value}`},[g("Deposited")],8,xs)),[[_]]),a("div",null,y(v.value||"-"),1)]),a("div",zs,[C((e(),t("div",{"data-copy":`Released: ${c.value}`},[g("Released")],8,Rs)),[[_]]),a("div",null,y(c.value||"-"),1)]),a("div",js,[C((e(),t("div",{"data-copy":`Authors: ${((B=n.value)==null?void 0:B.join(", "))||"-"}`},[g("Authors")],8,Ps)),[[_]]),n.value&&n.value.length?(e(),t("div",Ws,[(e(!0),t(D,null,F(n.value,(V,N)=>(e(),t(D,{key:N},[a("a",{href:z(V),target:"_blank",class:"lookup"},y(V),9,Us),n.value&&N<n.value.length-1?(e(),t("span",Ks,", ")):p("",!0)],64))),128))])):(e(),t("div",Gs,"-"))]),a("div",Js,[C((e(),t("div",{"data-copy":`Keywords: ${((q=r.value)==null?void 0:q.join(", "))||"-"}`},[g("Keywords")],8,qs)),[[_]]),r.value&&r.value.length?(e(),t("div",Qs,[(e(!0),t(D,null,F(r.value,(V,N)=>(e(),t(D,{key:N},[a("a",{href:z(V),target:"_blank",class:"lookup"},y(V),9,Xs),N<r.value.length-1?(e(),t("span",Ys,", ")):p("",!0)],64))),128))])):(e(),t("div",Zs,"-"))]),s.value?(e(),t("div",ea,[C((e(),t("div",{"data-copy":`PDB Entry: https://www.rcsb.org/structure/${u.value}`},[g("3D Data")],8,ta)),[[_]]),a("div",null,y($(d).mmolData3DFormat),1)])):p("",!0)]),oa,a("div",sa,[aa,a("div",la,[na,a("div",null,y(((X=(Q=o.value)==null?void 0:Q.Experimental)==null?void 0:X.Method)||"-"),1)]),a("div",ia,[ra,a("div",null,y(((Z=(Y=l.value)==null?void 0:Y.reflns_shell)==null?void 0:Z.d_res_high)||"-"),1)]),a("div",ca,[da,a("div",null,y(((te=(ee=l.value)==null?void 0:ee.refine)==null?void 0:te.ls_R_factor_R_free)||"-"),1)]),a("div",ua,[va,a("div",null,y(((se=(oe=l.value)==null?void 0:oe.refine)==null?void 0:se.ls_R_factor_R_work)||"-"),1)]),a("div",_a,[ma,a("div",null,y(((le=(ae=l.value)==null?void 0:ae.refine)==null?void 0:le.ls_R_factor_obs)||"-"),1)])]),ha,b.value?(e(),t("div",pa,[fa,a("img",{id:"validation-report",src:b.value,alt:"Validation report"},null,8,ya)])):p("",!0),b.value?(e(),t("hr",ga)):p("",!0),a("div",ka,[ba,(e(!0),t(D,null,F(o.value,(V,N)=>(e(),t(D,{key:N},[Object.keys(V).length?(e(),t("details",wa,[a("summary",null,[a("h4",null,[g(y(N),1),a("a",{href:"#"+$(d).mmolDataKeyMap[N],name:$(d).mmolDataKeyMap[N]},"#",8,$a)])]),a("div",null,[V.matrices&&V.vectors||Array.isArray(V)&&V[0].matrices&&V[0].vectors?(e(),M(ks,{key:0,data:V},null,8,["data"])):Array.isArray(V)?(e(),t("div",Sa,[T(ue,{data:V,allowCopy:!0},null,8,["data"])])):(e(),t("div",Ma,[T(ue,{data:V,allowCopy:!1},null,8,["data"])]))])])):p("",!0)],64))),128))])],64)}}}),Ca=H(Da,[["__scopeId","data-v-5c5b6fff"]]);const Ba={id:"content-wrap"},Ia={class:"col-left"};const Va=E({__name:"MolViewer",props:{context:{default:"file"},loading:{type:Boolean},loadingError:{}},emits:["retryLoad"],setup(h){const d=_e(),i=vt(),s=j(),l=J(),o=ce(),u=ye(),m=h,v=k(()=>s.molType),c=k(()=>{var r;return m.context=="molset"?l.breadCrumbPathArray.concat(["mol #"+((r=s.molFromMolsetIndex)==null?void 0:r.toString())]):l.breadCrumbPathArray});if(s.isEmpty&&l.data&&l.moduleName=="MolViewer"){const r=l.data;s.setMolData(r,"smol")}(m.context=="file"||m.context=="identifier")&&(s.inchi?s.fetchSmolVizData(s.inchi):s.smiles&&s.fetchSmolVizData(s.smiles)),$t(s.clear),window.onbeforeunload=function(){if(o.hasChanges)return!0;o.clear()},St(n),Mt(n);async function n(r,b,w){s.molFromMolset?w():s.hasChanges?await u.alert("If you leave, all changes to this molecule will be lost.",{title:"Unsaved molecule changes",primaryBtn:"Stay",secondaryBtn:"Discard",onCancel:()=>{r.fullPath!=b.fullPath&&s.clear(),w()},onSubmit:()=>w(!1)}):(r.path!=b.path&&s.clear(),w())}return(r,b)=>(e(),t(D,null,[a("div",{id:"mol-render",class:O({headless:$(i).headless})},[v.value=="smol"?(e(),M(Ft,{key:0})):v.value=="mmol"?(e(),M(Ht,{key:1})):p("",!0)],2),a("div",Ba,[a("div",Ia,[r.context=="identifier"?(e(),M(Vt,{key:0,backto:{name:"mol"}},{default:de(()=>[T(ne,{icon:"icn-file-json",iconHover:"icn-file-json-hover",btnStyle:"soft",mini:"",onClick:b[0]||(b[0]=w=>$(d).push("?use=json"))})]),_:1})):(e(),M(Dt,{key:1,pathArray:c.value},{default:de(()=>[T(ne,{icon:"icn-file-json",iconHover:"icn-file-json-hover",btnStyle:"soft",mini:"",onClick:b[1]||(b[1]=w=>$(d).push("?use=json"))})]),_:1},8,["pathArray"])),T(lo,{context:r.context,loading:r.loading},null,8,["context","loading"]),v.value=="smol"?(e(),M(rs,{key:2,loading:r.loading,loadingError:r.loadingError,onRetryLoad:b[2]||(b[2]=w=>r.$emit("retryLoad"))},null,8,["loading","loadingError"])):v.value=="mmol"?(e(),M(Ca,{key:3})):p("",!0)]),(r.loading,p("",!0))])],64))}}),Ta=H(Va,[["__scopeId","data-v-b7a671e5"]]),La=Object.freeze(Object.defineProperty({__proto__:null,default:Ta},Symbol.toStringTag,{value:"Module"}));export{Zt as B,Ta as M,Kt as a,Vt as b,La as c,mt as u};
