import{af as h}from"./index-DRw_S7Oy.js";function g(t,r){var e=Object.keys(t);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(t);r&&(n=n.filter(function(i){return Object.getOwnPropertyDescriptor(t,i).enumerable})),e.push.apply(e,n)}return e}function f(t){for(var r=1;r<arguments.length;r++){var e=arguments[r]!=null?arguments[r]:{};r%2?g(Object(e),!0).forEach(function(n){w(t,n,e[n])}):Object.getOwnPropertyDescriptors?Object.defineProperties(t,Object.getOwnPropertyDescriptors(e)):g(Object(e)).forEach(function(n){Object.defineProperty(t,n,Object.getOwnPropertyDescriptor(e,n))})}return t}function w(t,r,e){return r=x(r),r in t?Object.defineProperty(t,r,{value:e,enumerable:!0,configurable:!0,writable:!0}):t[r]=e,t}function y(t,r){if(t==null)return{};var e={},n=Object.keys(t),i,o;for(o=0;o<n.length;o++)i=n[o],!(r.indexOf(i)>=0)&&(e[i]=t[i]);return e}function b(t,r){if(t==null)return{};var e=y(t,r),n,i;if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(t);for(i=0;i<o.length;i++)n=o[i],!(r.indexOf(n)>=0)&&Object.prototype.propertyIsEnumerable.call(t,n)&&(e[n]=t[n])}return e}function O(t,r){if(typeof t!="object"||t===null)return t;var e=t[Symbol.toPrimitive];if(e!==void 0){var n=e.call(t,r||"default");if(typeof n!="object")return n;throw new TypeError("@@toPrimitive must return a primitive value.")}return(r==="string"?String:Number)(t)}function x(t){var r=O(t,"string");return typeof r=="symbol"?r:String(r)}var P=["width","height","viewBox"],j=["tabindex"],S={focusable:"false",preserveAspectRatio:"xMidYMid meet"};function _(){var t=arguments.length>0&&arguments[0]!==void 0?arguments[0]:{},r=t.width,e=t.height,n=t.viewBox,i=n===void 0?"0 0 ".concat(r," ").concat(e):n,o=b(t,P),l=o.tabindex,a=b(o,j),s=f(f(f({},S),a),{},{width:r,height:e,viewBox:i});return s["aria-label"]||s["aria-labelledby"]||s.title?(s.role="img",l!=null&&(s.focusable="true",s.tabindex=l)):s["aria-hidden"]=!0,s}function v(t){return h[t]}const p=v("h"),A=v("createApp"),d=(t,r,e)=>_({...r,preserveAspectRatio:"xMidYMid meet",xmlns:"http://www.w3.org/2000/svg",title:t,...e}),C=(t,r,e,n)=>{const i={attrs:d(t,r,e.attrs),on:n,style:{...e.staticStyle,...e.style}};return delete i.attrs.style,(e.staticClass||e.class)&&(i.class={},e.staticClass&&(i.class[e.staticClass]=!0),e.class&&(i.class[e.class]=!0)),i},B=(t,r,e)=>({props:{title:String},name:t,...A?{setup(n,i){let{title:o}=n,{attrs:l,slots:a}=i;return()=>p("svg",d(o,r,l),[...o?[p("title",o)]:[],...e.map(s=>{let{elem:u,attrs:c}=s;return p(u,c)}),...a.default?a.default():[]])}}:{functional:!0,render(n,i){let{props:{title:o},children:l,data:a,listeners:s}=i;return n("svg",C(o,r,a,s),[...o?[n("title",null,o)]:[],...e.map(u=>{let{elem:c,attrs:m}=u;return n(c,{attrs:m})}),...l||[]])}}}),D=B("Close16",{xmlns:"http://www.w3.org/2000/svg",viewBox:"0 0 32 32",fill:"currentColor",width:16,height:16},[{elem:"path",attrs:{d:"M17.4141 16L24 9.4141 22.5859 8 16 14.5859 9.4143 8 8 9.4141 14.5859 16 8 22.5859 9.4143 24 16 17.4141 22.5859 24 24 22.5859 17.4141 16z"}}]);export{D as C,B as c};
