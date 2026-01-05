import { Canvg } from 'https://cdn.skypack.dev/canvg@^4.0.0';
function flag_to_png(svg_string){
    const canvas = document.getElementById('png_flag');
    const ctx = canvas.getContext('2d');
    const v = Canvg.fromString(ctx, svg_string);
    v.start();
}
    import {data} from './usf_data.js'
        
          // jQuery.noConflict();
          // jQuery(function() {
          //   jQuery(".notetext").hide();
          //   jQuery(".note").click(function(event) {
          //     jQuery(this.nextSibling).toggle();
          //     event.stopPropagation();
          //   });
          //   jQuery("body").click(function(event) {
          //     jQuery(".notetext").hide();
          //   });
          // });
    function make_flag(x,y,method,alter,str_cou,canton_width,canton_height,aspect_ratio,first_stripe_red,star_points,star_points_2){
        let red_color = document.getElementById("red_color").value;
        let white_color = document.getElementById("white_color").value;
        let blue_color = document.getElementById("blue_color").value;
        let star_color = document.getElementById("star_color").value;
        var star_checked = null; 
        var inputElements = document.getElementById('star_color_different');
        if(inputElements.checked){
            star_checked = inputElements.value;
        }
        if(!star_checked){
            star_color = white_color;
        }
        let w = '';
        let adjusted_width = 19.0*aspect_ratio;
        let star_string = '<path id="s_pre" fill="#'+star_color+'" d="M';

        var inner_ring = Math.cos((star_points_2*Math.PI)/star_points)/Math.cos(Math.PI/star_points - (star_points_2*Math.PI)/star_points);
        for(let i=0;i<star_points;i++){
            var current_angle = Math.PI*2/star_points*i+Math.PI*3/2;
            var next_angle = Math.PI*2/star_points*(i+0.5)+Math.PI*3/2;
            star_string += (Math.cos(current_angle).toString()+','+Math.sin(current_angle).toString()+'L'+(inner_ring*Math.cos(next_angle)).toString()+','+(inner_ring*Math.sin(next_angle)).toString());
            if(i<star_points-1){
                star_string += 'L';
            }
        }
        star_string += ('z"/>');

        if(method == 0){
            let gap_limit = 0.58515/(19*7/13/12);
            let star_scale = gap_limit*Math.min(19*canton_height/y/2,adjusted_width*canton_width/x/2);
            let canton_height_svg = 19*canton_height;
            let canton_width_svg = adjusted_width*canton_width;


            w += ('<?xml version="1.0" encoding="UTF-8"?>\n<svg viewBox="0 0 '+(adjusted_width).toString()+' 19" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n<defs>\n'+star_string+'\n<g id="s"><use xlink:href="#s_pre" transform="scale('+star_scale.toString()+')"/></g>\n<g id="s'+x.toString()+'">');
            for(let i=0;i<x;i++){
                w += ('\n\t<use xlink:href="#s" x="'+(canton_width_svg/x/2*(2*i+1)).toString()+'" y="'+(19*canton_height/y/2).toString()+'"/>');
            }
            if(alter > 0){
                w += ('\n</g>\n\n<g id="s'+(x-1).toString()+'">');
                for(let i=0;i<x-1;i++){
                    w +=('\n\t<use xlink:href="#s" x="'+(canton_width_svg/x/2*(2*i+2)).toString()+'" y="'+(19*canton_height/y/2).toString()+'"/>');
                }
            }
            w += ('\n</g>\n<g id="u">\n\t<use xlink:href="#s'+(x).toString()+'"/>')
            for(let i=1;i<y;i++){
                if(alter > 0 && i >= (y-alter)/2 && i < (y+alter)/2){
                    w += ('\n\t<use xlink:href="#s'+(x-1).toString()+'" y="'+(19*canton_height/y/2*(2*i)).toString()+'"/>');
                }else{
                    w += ('\n\t<use xlink:href="#s'+(x).toString()+'" y="'+(19*canton_height/y/2*(2*i)).toString()+'"/>');
                }
            }
            w += ('</g>\n<rect id="stripe" width="'+(adjusted_width).toString()+'" height="'+(19/str_cou).toString()+'" y="'+((1-first_stripe_red)*19/str_cou).toString()+'" fill="#'+red_color+'"/>\n</defs>\n<rect width="'+(adjusted_width).toString()+'" height="19" fill="#'+white_color+'"/><use xlink:href="#stripe"/>');
            for(let i=2;i<str_cou-1+first_stripe_red;i+=2){
                w += ('<use xlink:href="#stripe" y="'+(i*19/str_cou).toString()+'"/>');
            }
            w += ('<rect width="'+(650*aspect_ratio*canton_width)/(650/19)+'" height="'+(650*canton_height)/(650/19)+'" fill="#'+blue_color+'"/><use xlink:href="#u"/></svg>');
            // w += ('\n</g>\n</defs>\n<rect width="36.1" height="19" fill="#'+red_color+'"/>\n<path stroke="#'+white_color+'" stroke-width="1.4615" d="   M0,2.19225H36.1   M0,5.11525H36.1   M0,8.03825H36.1   M0,10.96125H36.1   M0,13.88425H36.1   M0,16.80725H36.1"/>\n<rect width="14.44" height="10.23" fill="#'+blue_color+'"/>\n<use xlink:href="#u"/>\n</svg>');
        }else if(method == 1){
            let gap_limit = 100*0.0616/(7/13/5);
            let star_scale = gap_limit*Math.min(aspect_ratio*canton_width/(x+1)*2,canton_height/(y+1)*2)/2;

            w += ('<?xml version="1.0" encoding="UTF-8"?>\n<svg viewBox="0 0 '+(aspect_ratio*650).toString()+' 650" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n<defs>\n<polygon id="pt" points="-0.1624598481164531,0 0,-0.5 0.1624598481164531,0" transform="scale('+(star_scale).toString()+')" fill="#'+star_color+'"/>\n'+star_string+'\n<g id="s"><use xlink:href="#s_pre" transform="scale('+star_scale.toString()+')"/></g>\n<g id="s'+(Math.floor((x+1)/2)).toString()+'">');
            for(let i=0;i<Math.floor((x+1)/2);i++){
                w += ('\n\t<use xlink:href="#s" x="'+(100*aspect_ratio*canton_width/(x+1)*(2*i+1)).toString()+'"/>');
            }
            if(x%2 == 1){
                w += ('\n</g>\n<g id="s'+(Math.floor((x-1)/2)).toString()+'">');
                for(let i=0;i<Math.floor((x-1)/2);i++){
                    w += ('\n\t<use xlink:href="#s" x="'+(100*aspect_ratio*canton_width/(x+1)*(2*i+1)).toString()+'"/>');
                }
            }
            w += ('\n</g>\n<g id="u">');
            if(alter == 0){
                for(let i=0;i<Math.floor((y+1)/2);i++){
                    w += ('\n\t<use xlink:href="#s'+(Math.floor((x+1)/2)).toString()+'" y="'+(100*canton_height/(y+1)*(2*i+1)).toString()+'"/>');
                }
                for(let i=0;i<Math.floor((y)/2);i++){
                    w += ('\n\t<use xlink:href="#s'+(Math.floor((x)/2)).toString()+'" x="'+(100*aspect_ratio*canton_width/(x+1)*(1)).toString()+'" y="'+(100*canton_height/(y+1)*(2*i+2)).toString()+'"/>');
                }
            }else{
                let longer_row_count = Math.floor((y+1)/2)-alter;
                let shorter_row_count = Math.floor(y/2)+alter;

                var longer_rows_on_top = Math.floor((longer_row_count+1)/2);
                var longer_rows_on_bottom = Math.floor(longer_row_count/2);

                var row_types = [];
                for(let i=0;i<y;i++){
                    row_types.push(1-(i%2));
                }
                var longer_rows_current = Math.floor((y+1)/2);
                for(let i=0;i<Math.floor(y/2);i++){
                    if(longer_rows_current>longer_row_count){
                        if(row_types[Math.floor(y/2)-i]==1){
                            row_types[Math.floor(y/2)-i]=0;
                            longer_rows_current -= 1;
                        }
                    }
                    if(longer_rows_current>longer_row_count){
                        if(row_types[Math.floor(y/2)+i]==1){
                            row_types[Math.floor(y/2)+i]=0;
                            longer_rows_current -= 1;
                        }
                    }
                }
                for(let i=0;i<y;i++){
                    if(row_types[i]==0){
                        w += ('\n\t<use xlink:href="#s'+(Math.floor((x)/2)).toString()+'" x="'+(100*aspect_ratio*canton_width/(x+1)*(1)).toString()+'" y="'+(100*canton_height/(y+1)*(i+1)).toString()+'"/>');
                    }else{
                        w += ('\n\t<use xlink:href="#s'+(Math.floor((x+1)/2)).toString()+'" y="'+(100*canton_height/(y+1)*(i+1)).toString()+'"/>');
                    }
                }
            }
            w += ('</g>\n<rect id="stripe" width="'+(650*aspect_ratio).toString()+'" height="'+(650/str_cou).toString()+'" y="'+((1-first_stripe_red)*650/str_cou).toString()+'" fill="#'+red_color+'"/>\n</defs>\n<rect width="'+(650*aspect_ratio).toString()+'" height="650" fill="#'+white_color+'"/><use xlink:href="#stripe"/>');
            for(let i=2;i<str_cou-1+first_stripe_red;i+=2){
                w += ('<use xlink:href="#stripe" y="'+(i*650/str_cou).toString()+'"/>');
            }
            w += ('<rect width="'+(650*aspect_ratio*canton_width).toString()+'" height="'+(650*canton_height).toString()+'" fill="#'+blue_color+'"/><use xlink:href="#u" transform="scale(6.50)"/></svg>');
        }else if(method == 2){
            let gap_limit = 100*0.0616/(7/13/5);
            let star_scale = gap_limit*Math.min(aspect_ratio*canton_width/(x+1)*2,canton_height/(y+1)*2)/2;

            w += ('<?xml version="1.0" encoding="UTF-8"?>\n<svg viewBox="0 0 '+(aspect_ratio*650).toString()+' 650" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n<defs>\n<polygon id="pt" points="-0.1624598481164531,0 0,-0.5 0.1624598481164531,0" transform="scale('+(star_scale).toString()+')" fill="#'+star_color+'"/>\n'+star_string+'\n<g id="s"><use xlink:href="#s_pre" transform="scale('+star_scale.toString()+')"/></g>\n<g id="s'+(Math.floor((x+1)/2)).toString()+'">');
            for(let i=0;i<Math.floor((x+1)/2);i++){
                w += ('\n\t<use xlink:href="#s" x="'+(100*aspect_ratio*canton_width/(x+1)*(2*i+1)).toString()+'"/>');
            }
            if(x%2 == 1){
                w += ('\n</g>\n<g id="s'+(Math.floor((x-1)/2)).toString()+'">');
                for(let i=0;i<Math.floor((x-1)/2);i++){
                    w += ('\n\t<use xlink:href="#s" x="'+(100*aspect_ratio*canton_width/(x+1)*(2*i+1)).toString()+'"/>');
                }
            }
            w += ('\n</g>\n<g id="u">');
            for(let i=0;i<Math.floor((y+1)/2);i++){
                w += ('\n\t<use xlink:href="#s'+(Math.floor((x)/2)).toString()+'" x="'+(100*aspect_ratio*canton_width/(x+1)*(1)).toString()+'" y="'+(100*canton_height/(y+1)*(2*i+1)).toString()+'"/>');
            }
            for(let i=0;i<Math.floor((y)/2);i++){
                w += ('\n\t<use xlink:href="#s'+(Math.floor((x+1)/2)).toString()+'" y="'+(100*canton_height/(y+1)*(2*i+2)).toString()+'"/>');
            }
            w += ('</g>\n<rect id="stripe" width="'+(650*aspect_ratio).toString()+'" height="'+(650/str_cou).toString()+'" y="'+((1-first_stripe_red)*650/str_cou).toString()+'" fill="#'+red_color+'"/>\n</defs>\n<rect width="'+(650*aspect_ratio).toString()+'" height="650" fill="#'+white_color+'"/><use xlink:href="#stripe"/>');
            for(let i=2;i<str_cou-1+first_stripe_red;i+=2){
                w += ('<use xlink:href="#stripe" y="'+(i*650/str_cou).toString()+'"/>');
            }
            w += ('<rect width="'+(650*aspect_ratio*canton_width).toString()+'" height="'+(650*canton_height).toString()+'" fill="#'+blue_color+'"/><use xlink:href="#u" transform="scale(6.50)"/></svg>');
            // w += ('</g>\n<rect id="stripe" width="1235" height="50" fill="#'+red_color+'"/>\n</defs>\n<rect width="1235" height="650" fill="#'+white_color+'"/><use xlink:href="#stripe"/><use xlink:href="#stripe" y="100"/><use xlink:href="#stripe" y="200"/><use xlink:href="#stripe" y="300"/><use xlink:href="#stripe" y="400"/><use xlink:href="#stripe" y="500"/><use xlink:href="#stripe" y="600"/><rect width="494" height="350" fill="#'+blue_color+'"/><use xlink:href="#u" transform="scale(6.50)"/></svg>');
        }else if(method == 3){
            let gap_limit = 0.58515/(19*7/13/12);
            let canton_height_svg = 19*canton_height;
            let canton_width_svg = adjusted_width*canton_width;
            //let star_scale = gap_limit*3*Math.min(canton_width_svg/2,19*canton_height/2)/x;

            let star_scale = gap_limit*3*(Math.min(canton_width_svg/2,19*canton_height/2))/x/(1+gap_limit*3/0.8/x);
            

            w += ('<?xml version="1.0" encoding="UTF-8"?>\n<svg viewBox="0 0 '+(adjusted_width).toString()+' 19" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n<defs>\n'+star_string+'<g id="s"><use xlink:href="#s_pre" transform="scale('+star_scale.toString()+')"/></g>\n<g id="u">');
            for(let i=0;i<x;i++){
                w += ('\n\t<use xlink:href="#s" transform="translate('+(canton_width_svg/2).toString()+' '+(19*canton_height/2).toString()+'),rotate('+(i/x*360).toString()+'),translate(0 '+(-0.8*Math.min(canton_width_svg/2,19*canton_height/2)+star_scale).toString()+')"/>');
            }
            w += ('</g>')
            // if(alter > 0){
            //     w += ('\n</g>\n\n<g id="s'+(x-1).toString()+'">');
            //     for(let i=0;i<x-1;i++){
            //         w +=('\n\t<use xlink:href="#s" x="'+(canton_width_svg/x/2*(2*i+2)).toString()+'" y="'+(19*canton_height/y/2).toString()+'"/>');
            //     }
            // }
            // w += ('\n</g>\n<g id="u">\n\t<use xlink:href="#s'+(x).toString()+'"/>')
            // for(let i=1;i<y;i++){
            //     if(alter > 0 && i >= (y-alter)/2 && i < (y+alter)/2){
            //         w += ('\n\t<use xlink:href="#s'+(x-1).toString()+'" y="'+(19*canton_height/y/2*(2*i)).toString()+'"/>');
            //     }else{
            //         w += ('\n\t<use xlink:href="#s'+(x).toString()+'" y="'+(19*canton_height/y/2*(2*i)).toString()+'"/>');
            //     }
            // }
            w += ('\n<rect id="stripe" width="'+(adjusted_width).toString()+'" height="'+(19/str_cou).toString()+'" y="'+((1-first_stripe_red)*19/str_cou).toString()+'" fill="#'+red_color+'"/>\n</defs>\n<rect width="'+(adjusted_width).toString()+'" height="19" fill="#'+white_color+'"/><use xlink:href="#stripe"/>');
            for(let i=2;i<str_cou-1+first_stripe_red;i+=2){
                w += ('<use xlink:href="#stripe" y="'+(i*19/str_cou).toString()+'"/>');
            }
            w += ('<rect width="'+(650*aspect_ratio*canton_width)/(650/19)+'" height="'+(650*canton_height)/(650/19)+'" fill="#'+blue_color+'"/><use xlink:href="#u"/></svg>');
            // w += ('\n</g>\n</defs>\n<rect width="36.1" height="19" fill="#'+red_color+'"/>\n<path stroke="#'+white_color+'" stroke-width="1.4615" d="   M0,2.19225H36.1   M0,5.11525H36.1   M0,8.03825H36.1   M0,10.96125H36.1   M0,13.88425H36.1   M0,16.80725H36.1"/>\n<rect width="14.44" height="10.23" fill="#'+blue_color+'"/>\n<use xlink:href="#u"/>\n</svg>');
        }
        return w;
    }
    var svg_string = '';
    var star_count = 50;
    var stripe_count = 13;
    var aspect_ratio = parseFloat(document.getElementById("ratio_width").value)/parseFloat(document.getElementById("ratio_height").value);
    function modify_flag(){
        //document.write(data[parseInt(document.getElementById("star_count").value)][0]);
        star_count = parseInt(document.getElementById("star_count").value);
        stripe_count = parseInt(document.getElementById("stripe_count").value);
        var inputs = data[star_count].slice();
        var checkedValue = null; 
        var inputElements = document.getElementById('canton_height_units');
        if(inputElements.checked){
            checkedValue = inputElements.value;
        }
        var firstStripeRedValue = null; 
        var inputElements = document.getElementById('first_stripe_red');
        if(inputElements.checked){
            firstStripeRedValue = inputElements.value;
        }
        var inCircle = null; 
        var inputElements = document.getElementById('stars_in_circle');
        if(inputElements.checked){
            inCircle = inputElements.value;
        }

        let canton_height = 0;
        if (checkedValue){
            canton_height = parseFloat(document.getElementById("canton_height").value)/stripe_count;
        }else{
            canton_height = parseFloat(document.getElementById("canton_height").value)/100;

        }
        let canton_width = 0;
        canton_width = parseFloat(document.getElementById("canton_width").value)/100;

        aspect_ratio = parseFloat(document.getElementById("ratio_width").value)/parseFloat(document.getElementById("ratio_height").value);
        let star_points = 5;
        star_points = parseFloat(document.getElementById("star_points").value);
        let star_points_2 = 5;
        star_points_2 = parseFloat(document.getElementById("star_points_2").value);

        let first_stripe_red = 0;
        if (firstStripeRedValue){
            first_stripe_red = 1;
        }else{
            first_stripe_red = 0;
        }
        if(inCircle){
            inputs[0] = star_count;
            inputs[2] = 3;
        }
        svg_string = make_flag(inputs[0],inputs[1],inputs[2],inputs[3],stripe_count,canton_width,canton_height,aspect_ratio,first_stripe_red,star_points,star_points_2);
        var doc = new DOMParser().parseFromString(svg_string, 'application/xml');
        var el = document.getElementById("svg_flag");
        el.innerHTML = svg_string;
        //el.appendChild(
        //    el.ownerDocument.importNode(doc.documentElement, true)
        //);
    }
    function save(filename, data) {
        const blob = new Blob([data], {type: 'text/csv'});
        if(window.navigator.msSaveOrOpenBlob) {
            window.navigator.msSaveBlob(blob, filename);
        }
        else{
            const elem = window.document.createElement('a');
            elem.href = window.URL.createObjectURL(blob);
            elem.download = filename;        
            document.body.appendChild(elem);
            elem.click();        
            document.body.removeChild(elem);
        }
    }
    function download_flag(){
        save('us_'+(star_count).toString()+'_'+(stripe_count).toString()+'.svg', svg_string.replace('width="100%" height="52.63%"', 'width="1235" height="650"'));
    }
    function download_flag_png(){
        document.getElementById('canvas_hold').innerHTML = '<canvas id="png_flag" />';
        const canvas = document.getElementById('png_flag');
        canvas.height = Math.round(parseFloat(document.getElementById("png_height").value));
        canvas.width = Math.round(canvas.height*aspect_ratio);
        console.log(aspect_ratio)
        console.log(canvas.height)
        console.log(canvas.width)
        flag_to_png(svg_string);
        var anchor = document.createElement("a");
        anchor.href = canvas.toDataURL("image/png");
        anchor.download = 'us_'+(star_count).toString()+'_'+(stripe_count).toString()+'.png';
        anchor.click();
        document.getElementById('canvas_hold').innerHTML = '';
        // canvas.toBlob(function(blob) {
        //     saveAs(blob, 'us_'+(star_count).toString()+'_'+(stripe_count).toString()+'.png');
        // }, "image/png");
        //save('us_'+(star_count).toString()+'_'+(stripe_count).toString()+'.png', );
    }

    

    export {modify_flag, download_flag, download_flag_png};