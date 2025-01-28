import{d as w,L as g,h as _,i as C,v as e,j as i,q as l,n as t,I as o,Q as s,F as k,G as h,C as y,D as S,E as $}from"./index-DwYRJoEC.js";import{c as B,C as H}from"./16-Oa6gYw0C.js";const I=B("ChevronRight16",{xmlns:"http://www.w3.org/2000/svg",viewBox:"0 0 16 16",fill:"currentColor",width:16,height:16},[{elem:"path",attrs:{d:"M11 8L6 13 5.3 12.3 9.6 8 5.3 3.7 6 3z"}}]),c=r=>(y("data-v-b0760941"),r=r(),S(),r),M=c(()=>e("hr",null,null,-1)),T=c(()=>e("h4",null,":: Modals",-1)),x=c(()=>e("br",null,null,-1)),F=c(()=>e("br",null,null,-1)),N=h("<br data-v-b0760941><br data-v-b0760941><br data-v-b0760941><br data-v-b0760941><hr data-v-b0760941><h4 data-v-b0760941>:: Headers</h4><h1 data-v-b0760941>Header 1</h1><h2 data-v-b0760941>Header 2</h2><h3 data-v-b0760941>Header 3</h3><h4 data-v-b0760941>Header 4</h4><h5 data-v-b0760941>Header 5</h5><br data-v-b0760941><br data-v-b0760941><hr data-v-b0760941><h4 data-v-b0760941>:: Icons</h4>",15),V={style:{color:"red"}},z={href:"#"},E=c(()=>e("br",null,null,-1)),A=c(()=>e("br",null,null,-1)),L=c(()=>e("hr",null,null,-1)),j=c(()=>e("h4",null,":: Icons",-1)),O={class:"icons-wrap"},D={style:{background:"#ffd"}},G=c(()=>e("br",null,null,-1)),K={class:"icons-in-field"},P={class:"icons-wrap"},R=h('<span data-v-b0760941></span><br data-v-b0760941><br data-v-b0760941><hr data-v-b0760941><h4 data-v-b0760941>:: Colors</h4><div class="swatch-wrap" data-v-b0760941><div class="swatch black" data-v-b0760941>$black</div><div class="swatch black-60" data-v-b0760941>$black-60</div><div class="swatch black-30" data-v-b0760941>$black-30</div><div class="swatch black-20" data-v-b0760941>$black-20</div><div class="swatch black-10" data-v-b0760941>$black-10</div><div class="swatch soft-bg" data-v-b0760941>$soft-bg</div></div><div class="swatch-wrap" data-v-b0760941><div class="swatch blue" data-v-b0760941>$blue</div><div class="swatch blue-hover" data-v-b0760941>$blue-hover</div><div class="swatch blue-30" data-v-b0760941>$blue-30</div><div class="swatch blue-20" data-v-b0760941>$blue-20</div><div class="swatch blue-10" data-v-b0760941>$blue-10</div><div class="swatch blue-05" data-v-b0760941>$blue-05</div></div><div class="swatch-wrap" data-v-b0760941><div class="swatch success" data-v-b0760941>$success</div><div class="swatch warning" data-v-b0760941>$warning</div><div class="swatch caution" data-v-b0760941>$caution</div><div class="swatch error" data-v-b0760941>$error</div></div>',8),U=w({__name:"KitchenSink",setup(r){const n=g();async function f(){const u=await n.alert("Hello world",{title:"I block the thread",secondaryBtn:!0,onSubmit:v,onCancel:b});console.log(u?"Continue after SUBMIT":"Continue after CANCEL")}async function m(){const u=await n.confirm("Are you sure?",{onSubmit:v,onCancel:b});console.log(u?"Continue after SUBMIT":"Continue after CANCEL")}function v(){alert("yes"),n.hide()}function b(){alert("no"),n.hide()}function p(){alert("Other button"),n.hide()}return(u,a)=>(_(),C(k,null,[M,T,e("button",{onClick:a[0]||(a[0]=d=>i(n).alert("Hello world"))},"simple alert"),l("   "),e("button",{onClick:a[1]||(a[1]=d=>i(n).alert("Hello <span style='color: red'>world</span>",{html:!0,size:"md",primaryBtn:"One",secondaryBtn:"Two",otherBtn:"Three",title:"Here's a title",onSubmit:v,onCancel:b,onOther:p}))}," advanced alert "),l("   "),e("button",{onClick:f},"alert promise"),l("   "),e("button",{onClick:m},"confirm promise"),l("   "),e("button",{onClick:a[2]||(a[2]=d=>i(n).confirm("Are you sure?",{onSubmit:v,onCancel:b}))},"confirm modal"),x,F,e("button",{onClick:a[3]||(a[3]=d=>i(n).display("_ModalTemplate"))},"template"),l("  "),e("button",{onClick:a[4]||(a[4]=d=>i(n).display("ModalViewer"))},"file type"),l("  "),e("button",{onClick:a[5]||(a[5]=d=>i(n).display("ModalWorkspaces"))},"workspaces"),l("  "),e("button",{onClick:a[6]||(a[6]=d=>i(n).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",ext:"json",ext2:"molset"}))}," Save file "),l("   "),e("button",{onClick:a[7]||(a[7]=d=>i(n).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",dataType:"molset"}))}," Save molset as "),l("   "),e("button",{onClick:a[8]||(a[8]=d=>i(n).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",dataType:"smol"},{onSubmit:()=>{console.log("submitted")},onCancel:()=>{console.log("cancelled")}}))}," Save mol as "),l("   "),N,e("span",V,[e("a",z,[t(o,{icon:"icn-file-smol"})]),t(o,{icon:"icn-file-smol"}),t(o,{icon:"icn-file-mmol"}),t(o,{icon:"icn-file-molset"}),t(o,{icon:"icn-file-data"}),t(o,{icon:"icn-file-text"}),t(o,{icon:"icn-link"}),t(o,{icon:"icn-reaction"}),t(o,{icon:"icn-star-full"}),t(o,{icon:"icn-file-run"}),t(o,{icon:"icn-file-json"}),t(o,{icon:"icn-file-pdf"}),t(o,{icon:"icn-file-run"}),t(o,{icon:"icn-file-sdf"}),t(o,{icon:"icn-file-molset-csv"}),t(o,{icon:"icn-file-svg"}),t(o,{icon:"icn-file-md"}),t(o,{icon:"icn-file-unk"}),t(o,{icon:"icn-caret-left"}),t(o,{icon:"icn-caret-right"}),t(o,{icon:"icn-caret-up"}),t(o,{icon:"icn-caret-down"}),t(o,{icon:"icn-file-smol",size:"large"}),t(i(H)),t(i(I))]),E,A,L,j,e("div",O,[e("div",D,[l(" Default "),t(s,{icon:"icn-star-full"})]),e("div",null,[l(" Soft "),t(s,{icon:"icn-star-full",btnStyle:"soft"})]),e("div",null,[l(" Carbon "),t(s,{icon:"icn-star-full",btnStyle:"carbon"})]),e("div",null,[l(" Custom colors "),t(s,{icon:"icn-star-full",color:"green",colorHover:"red"})]),e("div",null,[l(" Toggle "),t(s,{icon:"icn-star-full",toggle:!0})]),e("div",null,[l(" Toggle with custom color "),t(s,{icon:"icn-star-full",toggle:!0,colorToggle:"#d3bf0b"})]),e("div",null,[l(" Mini "),t(s,{icon:"icn-star-full",mini:!0})])]),G,e("div",K,[l(" Examples "),e("div",P,[t(s,{icon:"icn-full-screen-large",iconHover:"icn-full-screen-large-hover",btnStyle:"soft",icnSize:"large"}),t(s,{icon:"icn-star-full-large",iconHover:"icn-star-full",iconSel:"icn-star-full",colorToggle:"#d3bf0b",toggle:!0,icnSize:"large"})])]),R],64))}}),W=$(U,[["__scopeId","data-v-b0760941"]]);export{W as default};
