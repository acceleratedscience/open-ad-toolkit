import{d as w,L as g,h as _,i as C,v as t,j as l,q as n,n as e,I as o,Q as s,F as k,G as h,C as S,D as y,E as $}from"./index-CbxKm7xL.js";import{c as B,C as H}from"./16-Clj_678a.js";const M=B("ChevronRight16",{xmlns:"http://www.w3.org/2000/svg",viewBox:"0 0 16 16",fill:"currentColor",width:16,height:16},[{elem:"path",attrs:{d:"M11 8L6 13 5.3 12.3 9.6 8 5.3 3.7 6 3z"}}]),c=r=>(S("data-v-28e7782e"),r=r(),y(),r),T=c(()=>t("hr",null,null,-1)),I=c(()=>t("h4",null,":: Modals",-1)),x=c(()=>t("br",null,null,-1)),F=c(()=>t("br",null,null,-1)),N=h("<br data-v-28e7782e><br data-v-28e7782e><br data-v-28e7782e><br data-v-28e7782e><hr data-v-28e7782e><h4 data-v-28e7782e>:: Headers</h4><h1 data-v-28e7782e>Header 1</h1><h2 data-v-28e7782e>Header 2</h2><h3 data-v-28e7782e>Header 3</h3><h4 data-v-28e7782e>Header 4</h4><h5 data-v-28e7782e>Header 5</h5><br data-v-28e7782e><br data-v-28e7782e><hr data-v-28e7782e><h4 data-v-28e7782e>:: Icons</h4>",15),V={style:{color:"red"}},z={href:"#"},E=c(()=>t("br",null,null,-1)),A=c(()=>t("br",null,null,-1)),L=c(()=>t("hr",null,null,-1)),j=c(()=>t("h4",null,":: Icons",-1)),O={class:"icons-wrap"},D={style:{background:"#ffd"}},G=c(()=>t("br",null,null,-1)),K={class:"icons-in-field"},P={class:"icons-wrap"},R=h('<span data-v-28e7782e></span><br data-v-28e7782e><br data-v-28e7782e><hr data-v-28e7782e><h4 data-v-28e7782e>:: Colors</h4><div class="swatch-wrap" data-v-28e7782e><div class="swatch black" data-v-28e7782e>$black</div><div class="swatch black-60" data-v-28e7782e>$black-60</div><div class="swatch black-30" data-v-28e7782e>$black-30</div><div class="swatch black-20" data-v-28e7782e>$black-20</div><div class="swatch black-10" data-v-28e7782e>$black-10</div><div class="swatch soft-bg" data-v-28e7782e>$soft-bg</div></div><div class="swatch-wrap" data-v-28e7782e><div class="swatch blue" data-v-28e7782e>$blue</div><div class="swatch blue-hover" data-v-28e7782e>$blue-hover</div><div class="swatch blue-30" data-v-28e7782e>$blue-30</div><div class="swatch blue-20" data-v-28e7782e>$blue-20</div><div class="swatch blue-10" data-v-28e7782e>$blue-10</div><div class="swatch blue-05" data-v-28e7782e>$blue-05</div></div><div class="swatch-wrap" data-v-28e7782e><div class="swatch success" data-v-28e7782e>$success</div><div class="swatch warning" data-v-28e7782e>$warning</div><div class="swatch caution" data-v-28e7782e>$caution</div><div class="swatch error" data-v-28e7782e>$error</div></div>',8),U=w({__name:"KitchenSink",setup(r){const i=g();async function f(){const b=await i.alert("Hello world",{title:"I block the thread",secondaryBtn:!0,onSubmit:v,onCancel:u});console.log(b?"Continue after SUBMIT":"Continue after CANCEL")}async function m(){const b=await i.confirm("Are you sure?",{onSubmit:v,onCancel:u});console.log(b?"Continue after SUBMIT":"Continue after CANCEL")}function v(){alert("yes"),i.hide()}function u(){alert("no"),i.hide()}function p(){alert("Other button"),i.hide()}return(b,a)=>(_(),C(k,null,[T,I,t("button",{onClick:a[0]||(a[0]=d=>l(i).alert("Hello world"))},"simple alert"),n("   "),t("button",{onClick:a[1]||(a[1]=d=>l(i).alert("Hello <span style='color: red'>world</span>",{html:!0,size:"md",primaryBtn:"One",secondaryBtn:"Two",otherBtn:"Three",title:"Here's a title",onSubmit:v,onCancel:u,onOther:p}))}," advanced alert "),n("   "),t("button",{onClick:f},"alert promise"),n("   "),t("button",{onClick:m},"confirm promise"),n("   "),t("button",{onClick:a[2]||(a[2]=d=>l(i).confirm("Are you sure?",{onSubmit:v,onCancel:u}))},"confirm modal"),x,F,t("button",{onClick:a[3]||(a[3]=d=>l(i).display("_ModalTemplate"))},"template"),n("  "),t("button",{onClick:a[4]||(a[4]=d=>l(i).display("ModalViewer"))},"file type"),n("  "),t("button",{onClick:a[5]||(a[5]=d=>l(i).display("ModalWorkspaces"))},"workspaces"),n("  "),t("button",{onClick:a[6]||(a[6]=d=>l(i).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",ext:"json",ext2:"molset"}))}," Save file "),n("   "),t("button",{onClick:a[7]||(a[7]=d=>l(i).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",dataType:"molset"}))}," Save molset as "),n("   "),t("button",{onClick:a[8]||(a[8]=d=>l(i).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",dataType:"smol"},{onSubmit:()=>{console.log("submitted")},onCancel:()=>{console.log("cancelled")}}))}," Save mol as "),n("   "),N,t("span",V,[t("a",z,[e(o,{icon:"icn-file-smol"})]),e(o,{icon:"icn-file-smol"}),e(o,{icon:"icn-file-mmol"}),e(o,{icon:"icn-file-molset"}),e(o,{icon:"icn-file-data"}),e(o,{icon:"icn-file-text"}),e(o,{icon:"icn-link"}),e(o,{icon:"icn-reaction"}),e(o,{icon:"icn-star"}),e(o,{icon:"icn-file-run"}),e(o,{icon:"icn-file-json"}),e(o,{icon:"icn-file-pdf"}),e(o,{icon:"icn-file-run"}),e(o,{icon:"icn-file-sdf"}),e(o,{icon:"icn-file-molset-csv"}),e(o,{icon:"icn-file-svg"}),e(o,{icon:"icn-file-md"}),e(o,{icon:"icn-file-unk"}),e(o,{icon:"icn-caret-left"}),e(o,{icon:"icn-caret-right"}),e(o,{icon:"icn-caret-up"}),e(o,{icon:"icn-caret-down"}),e(o,{icon:"icn-file-smol",size:"large"}),e(l(H)),e(l(M))]),E,A,L,j,t("div",O,[t("div",D,[n(" Default "),e(s,{icon:"icn-star"})]),t("div",null,[n(" Soft "),e(s,{icon:"icn-star",btnStyle:"soft"})]),t("div",null,[n(" Carbon "),e(s,{icon:"icn-star",btnStyle:"carbon"})]),t("div",null,[n(" Custom colors "),e(s,{icon:"icn-star",color:"green",colorHover:"red"})]),t("div",null,[n(" Toggle "),e(s,{icon:"icn-star",toggle:!0})]),t("div",null,[n(" Toggle with custom color "),e(s,{icon:"icn-star",toggle:!0,colorToggle:"#d3bf0b"})]),t("div",null,[n(" Mini "),e(s,{icon:"icn-star",mini:!0})])]),G,t("div",K,[n(" Examples "),t("div",P,[e(s,{icon:"icn-full-screen-large",iconHover:"icn-full-screen-large-hover",btnStyle:"soft",icnSize:"large"}),e(s,{icon:"icn-star-large-outline",iconHover:"icn-star",iconSel:"icn-star",colorToggle:"#d3bf0b",toggle:!0,icnSize:"large"})])]),R],64))}}),W=$(U,[["__scopeId","data-v-28e7782e"]]);export{W as default};