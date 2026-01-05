import { Canvg } from 'https://cdn.skypack.dev/canvg@^4.0.0';
function flag_to_png(svg_string){
	const canvas = document.getElementById('png_flag');
	const ctx = canvas.getContext('2d');
	const v = Canvg.fromString(ctx, svg_string);

}
export {flag_to_png};