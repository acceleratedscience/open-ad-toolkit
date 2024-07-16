import{d as w,K as g,i as _,j as k,x as e,k as l,v as n,p as t,H as a,I as c,F as C,P as h,C as S,D as y,E as $}from"./index-juGIK60X.js";import{c as H,C as B}from"./16-e0gwu24e.js";const M=H("ChevronRight16",{xmlns:"http://www.w3.org/2000/svg",viewBox:"0 0 16 16",fill:"currentColor",width:16,height:16},[{elem:"path",attrs:{d:"M11 8L6 13 5.3 12.3 9.6 8 5.3 3.7 6 3z"}}]),s=r=>(S("data-v-593c4083"),r=r(),y(),r),T=s(()=>e("hr",null,null,-1)),x=s(()=>e("h4",null,":: Modals",-1)),I=s(()=>e("br",null,null,-1)),F=s(()=>e("br",null,null,-1)),N=h("<br data-v-593c4083><br data-v-593c4083><br data-v-593c4083><br data-v-593c4083><hr data-v-593c4083><h4 data-v-593c4083>:: Headers</h4><h1 data-v-593c4083>Header 1</h1><h2 data-v-593c4083>Header 2</h2><h3 data-v-593c4083>Header 3</h3><h4 data-v-593c4083>Header 4</h4><h5 data-v-593c4083>Header 5</h5><br data-v-593c4083><br data-v-593c4083><hr data-v-593c4083><h4 data-v-593c4083>:: Icons</h4>",15),V={style:{color:"red"}},z={href:"#"},E=s(()=>e("br",null,null,-1)),A=s(()=>e("br",null,null,-1)),j=s(()=>e("hr",null,null,-1)),K=s(()=>e("h4",null,":: Icons",-1)),L={class:"icons-wrap"},O={style:{background:"#ffd"}},P=s(()=>e("br",null,null,-1)),D={class:"icons-in-field"},R={class:"icons-wrap"},U=h('<span data-v-593c4083></span><br data-v-593c4083><br data-v-593c4083><hr data-v-593c4083><h4 data-v-593c4083>:: Colors</h4><div class="swatch-wrap" data-v-593c4083><div class="swatch black" data-v-593c4083>$black</div><div class="swatch black-60" data-v-593c4083>$black-60</div><div class="swatch black-30" data-v-593c4083>$black-30</div><div class="swatch black-20" data-v-593c4083>$black-20</div><div class="swatch black-10" data-v-593c4083>$black-10</div><div class="swatch soft-bg" data-v-593c4083>$soft-bg</div></div><div class="swatch-wrap" data-v-593c4083><div class="swatch blue" data-v-593c4083>$blue</div><div class="swatch blue-hover" data-v-593c4083>$blue-hover</div><div class="swatch blue-30" data-v-593c4083>$blue-30</div><div class="swatch blue-20" data-v-593c4083>$blue-20</div><div class="swatch blue-10" data-v-593c4083>$blue-10</div><div class="swatch blue-05" data-v-593c4083>$blue-05</div></div><div class="swatch-wrap" data-v-593c4083><div class="swatch success" data-v-593c4083>$success</div><div class="swatch warning" data-v-593c4083>$warning</div><div class="swatch caution" data-v-593c4083>$caution</div><div class="swatch error" data-v-593c4083>$error</div></div>',8),G=w({__name:"KitchenSink",setup(r){const i=g();async function f(){const b=await i.alert("Hello world",{title:"I block the thread",secondaryBtn:!0,onSubmit:v,onCancel:u});console.log(b?"Continue after SUBMIT":"Continue after CANCEL")}async function m(){const b=await i.confirm("Are you sure?",{onSubmit:v,onCancel:u});console.log(b?"Continue after SUBMIT":"Continue after CANCEL")}function v(){alert("yes"),i.hide()}function u(){alert("no"),i.hide()}function p(){alert("Other button"),i.hide()}return(b,o)=>(_(),k(C,null,[T,x,e("button",{onClick:o[0]||(o[0]=d=>l(i).alert("Hello world"))},"simple alert"),n("   "),e("button",{onClick:o[1]||(o[1]=d=>l(i).alert("Hello <span style='color: red'>world</span>",{html:!0,size:"md",primaryBtn:"One",secondaryBtn:"Two",otherBtn:"Three",title:"Here's a title",onSubmit:v,onCancel:u,onOther:p}))}," advanced alert "),n("   "),e("button",{onClick:f},"alert promise"),n("   "),e("button",{onClick:m},"confirm promise"),n("   "),e("button",{onClick:o[2]||(o[2]=d=>l(i).confirm("Are you sure?",{onSubmit:v,onCancel:u}))},"confirm modal"),I,F,e("button",{onClick:o[3]||(o[3]=d=>l(i).display("_ModalTemplate"))},"template"),n("  "),e("button",{onClick:o[4]||(o[4]=d=>l(i).display("ModalViewer"))},"file type"),n("  "),e("button",{onClick:o[5]||(o[5]=d=>l(i).display("ModalWorkspaces"))},"workspaces"),n("  "),e("button",{onClick:o[6]||(o[6]=d=>l(i).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",ext:"json",ext2:"molset"}))}," Save file "),n("   "),e("button",{onClick:o[7]||(o[7]=d=>l(i).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",dataType:"molset"}))}," Save molset as "),n("   "),e("button",{onClick:o[8]||(o[8]=d=>l(i).display("ModalSaveFile",{path:"my_dir/sub_dir/subsub_dir",filename:"Foobar",dataType:"mol"},{onSubmit:()=>{console.log("submitted")},onCancel:()=>{console.log("cancelled")}}))}," Save mol as "),n("   "),N,e("span",V,[e("a",z,[t(a,{icon:"icn-file-mol"})]),t(a,{icon:"icn-file-mol"}),t(a,{icon:"icn-file-molset"}),t(a,{icon:"icn-file-data"}),t(a,{icon:"icn-file-text"}),t(a,{icon:"icn-link"}),t(a,{icon:"icn-reaction"}),t(a,{icon:"icn-star"}),t(a,{icon:"icn-file-run"}),t(a,{icon:"icn-file-json"}),t(a,{icon:"icn-file-pdf"}),t(a,{icon:"icn-file-run"}),t(a,{icon:"icn-file-sdf"}),t(a,{icon:"icn-file-molset-csv"}),t(a,{icon:"icn-file-svg"}),t(a,{icon:"icn-file-md"}),t(a,{icon:"icn-file-unk"}),t(a,{icon:"icn-caret-left"}),t(a,{icon:"icn-caret-right"}),t(a,{icon:"icn-caret-up"}),t(a,{icon:"icn-caret-down"}),t(a,{icon:"icn-file-mol",size:"large"}),t(l(B)),t(l(M))]),E,A,j,K,e("div",L,[e("div",O,[n(" Default "),t(c,{icon:"icn-star"})]),e("div",null,[n(" Soft "),t(c,{icon:"icn-star",btnStyle:"soft"})]),e("div",null,[n(" Carbon "),t(c,{icon:"icn-star",btnStyle:"carbon"})]),e("div",null,[n(" Custom colors "),t(c,{icon:"icn-star",color:"green",colorHover:"red"})]),e("div",null,[n(" Toggle "),t(c,{icon:"icn-star",toggle:!0})]),e("div",null,[n(" Toggle with custom color "),t(c,{icon:"icn-star",toggle:!0,colorToggle:"#d3bf0b"})]),e("div",null,[n(" Mini "),t(c,{icon:"icn-star",mini:!0})])]),P,e("div",D,[n(" Examples "),e("div",R,[t(c,{icon:"icn-full-screen-large",iconHover:"icn-full-screen-large-hover",btnStyle:"soft",icnSize:"large"}),t(c,{icon:"icn-star-large-outline",iconHover:"icn-star",iconSel:"icn-star",colorToggle:"#d3bf0b",toggle:!0,icnSize:"large"})])]),U],64))}}),J=$(G,[["__scopeId","data-v-593c4083"]]);export{J as default};