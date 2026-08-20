let FLAGS, FLAG_ORDER, CONSTELLATIONS, BAND_COLORS, STAR_DATA, CONST_LINES, CONST_LINE_RANGES;

let flag, rotmat, zoom, rollRad, pxPerUnit, centerX, centerY, clipRegion,
    flagW, flagH, reflect, magBounds, bands, magMinBase, magMaxBase, magMin, magMax,
    liveStarsGroup, starsPickedPath, constLinesPath,
    centerCrosshairG, crosshairOutline, crosshairFg;

let visibleStars = [];
let hiddenStars = new Set;
let hoveredStar = null;
let pickMode = false;
let moveCenterMode = false;
let isMovingCenter = false;
let showHiddenMode = false;
let showConstLines = false;
let onlyLinesMode = false;
let mirrorView = false;

const flagSvg = document.getElementById("flag");

async function flag_to_png(svg_string, pngW, pngH){
    const canvas = document.getElementById('png_flag');
    canvas.width = pngW;
    canvas.height = pngH;
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, pngW, pngH);
    const url = URL.createObjectURL(new Blob([svg_string], {type: "image/svg+xml;charset=utf-8"}));
    try {
        const img = new Image();
        await new Promise((resolve, reject) => {
            img.onload = resolve;
            img.onerror = reject;
            img.src = url;
        });
        ctx.drawImage(img, 0, 0, pngW, pngH);
    } finally {
        URL.revokeObjectURL(url);
    }
    return new Promise(resolve => canvas.toBlob(resolve, "image/png"));
}

function dot(a, b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

function cross(a, b){
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]];
}

function normalize(v){
    const m = Math.hypot(v[0], v[1], v[2]);
    return [v[0]/m, v[1]/m, v[2]/m];
}

function rotateVec(v, axis, angle){
    const c = Math.cos(angle), s = Math.sin(angle), d = dot(v, axis), k = cross(axis, v);
    return [
        v[0]*c + k[0]*s + axis[0]*d*(1 - c),
        v[1]*c + k[1]*s + axis[1]*d*(1 - c),
        v[2]*c + k[2]*s + axis[2]*d*(1 - c)
    ];
}

function rotateBasis(axis, angle){
    rotmat = rotmat.map(v => rotateVec(v, axis, angle));
}

function reorthonormalizeBasis(){
    const fwd = rotmat[2];
    let right = cross([0, 0, 1], fwd);
    if(Math.hypot(right[0], right[1], right[2]) < 1e-6) right = cross([0, 1, 0], fwd);
    right = normalize(right);
    let up = cross(fwd, right);
    if(dot(right, rotmat[0]) < 0) right = [-right[0], -right[1], -right[2]];
    if(dot(up, rotmat[1]) < 0) up = [-up[0], -up[1], -up[2]];
    rotmat = [right, up, fwd];
}

function rolledAxes(){
    const c = Math.cos(rollRad), s = Math.sin(rollRad);
    return {
        rRight: [0, 1, 2].map(i => c*rotmat[0][i] + s*rotmat[1][i]),
        rUp: [0, 1, 2].map(i => -s*rotmat[0][i] + c*rotmat[1][i]),
        fwd: rotmat[2]
    }
}

function starPath(cx, cy, outerR, innerR, points, rotationDeg){
    const rot = (rotationDeg || 0)*Math.PI/180;
    let d = "";
    for(let i = 0; i < points; i += 1){
        const outerA = -Math.PI/2 + rot + 2*i*Math.PI/points,
              innerA = -Math.PI/2 + rot + 2*(i + .5)*Math.PI/points,
              ox = cx + outerR*Math.cos(outerA), oy = cy + outerR*Math.sin(outerA),
              ix = cx + innerR*Math.cos(innerA), iy = cy + innerR*Math.sin(innerA);
        d += (i === 0 ? "M" : "L") + ox.toFixed(1) + "," + oy.toFixed(1) + "L" + ix.toFixed(1) + "," + iy.toFixed(1);
    }
    return d + "z";
}

function ratioOfSkip(points, skip){
    const a = Math.PI/points;
    return Math.cos(skip*a)/Math.cos((skip - 1)*a);
}

function skipOfRatio(points, ratio){
    const a = Math.PI/points;
    return Math.atan2(1 - ratio*Math.cos(a), ratio*Math.sin(a))/a;
}

function wrapSkip(skip, points){
    skip = (skip%points + points)%points;
    return skip > points/2 ? points - skip : skip;
}

function bandIndexForMag(mag){
    for(let i = 0; i < bands.length; i += 1)
        if(mag >= magBounds[i] && mag < magBounds[i + 1]) return i;
    return mag < magBounds[0] ? 0 : bands.length - 1;
}

let visibleMags = [];

const LINE_COS_THRESHOLD = Math.cos(0.3*Math.PI/180);

function constellationEndpoints(abbr){
    const range = CONST_LINE_RANGES[abbr];
    if(!range) return null;
    const [lo, hi] = range, pts = [];
    for(let i = lo; i < hi; i += 1){
        const o = 6*i;
        pts.push([CONST_LINES[o], CONST_LINES[o + 1], CONST_LINES[o + 2]]);
        pts.push([CONST_LINES[o + 3], CONST_LINES[o + 4], CONST_LINES[o + 5]]);
    }
    return pts;
}

function isOnConstellationLines(x, y, z, pts){
    if(!pts) return false;
    for(const p of pts){
        if(x*p[0] + y*p[1] + z*p[2] > LINE_COS_THRESHOLD) return true;
    }
    return false;
}

function render(){
    const scale = pxPerUnit*zoom, axes = rolledAxes(), fwd = axes.fwd, rUp = axes.rUp,
          rRight = mirrorView ? axes.rRight.map(v => -v) : axes.rRight;
    syncPointingFields(fwd);
    visibleMags = [];
    visibleStars = [];
    const starCount = STAR_DATA.length/4,
          lineEndpoints = onlyLinesMode ? constellationEndpoints(constellationSelect.value) : null;
    for(let i = 0; i < starCount; i += 1){
        const o = 4*i, x = STAR_DATA[o], y = STAR_DATA[o + 1], z = STAR_DATA[o + 2], mag = STAR_DATA[o + 3];
        if(fwd[0]*x + fwd[1]*y + fwd[2]*z <= 0) continue;
        if(onlyLinesMode && !isOnConstellationLines(x, y, z, lineEndpoints)) continue;
        const sx = centerX + scale*(rRight[0]*x + rRight[1]*y + rRight[2]*z),
              sy = centerY + scale*(rUp[0]*x + rUp[1]*y + rUp[2]*z),
              band = bands[bandIndexForMag(mag)], pad = band.size;
        if(sx < clipRegion[0] - pad || sx > clipRegion[2] + pad || sy < clipRegion[1] - pad || sy > clipRegion[3] + pad) continue;
        visibleMags.push(mag);
        if(band.render) visibleStars.push({j: i, sx, sy, b: band})
    }
    renderStarPaths();
    rebuildBarDots();
    renderConstLines(scale, rRight, rUp, fwd);
    updateCrosshair();
}

function updateCrosshair(){
    if(!crosshairFg) return;
    const arm = 20*flagW/760,
          d = `M${(centerX - arm).toFixed(2)},${centerY.toFixed(2)}L${(centerX + arm).toFixed(2)},${centerY.toFixed(2)}` +
              `M${centerX.toFixed(2)},${(centerY - arm).toFixed(2)}L${centerX.toFixed(2)},${(centerY + arm).toFixed(2)}`;
    crosshairOutline.setAttribute("d", d);
    crosshairFg.setAttribute("d", d);
}

function renderConstLines(scale, rRight, rUp, fwd){
    if(!constLinesPath) return;
    if(!showConstLines) return void constLinesPath.setAttribute("d", "");
    let d = "";
    const segCount = CONST_LINES.length/6;
    for(let i = 0; i < segCount; i += 1){
        const o = 6*i,
              x1 = CONST_LINES[o], y1 = CONST_LINES[o + 1], z1 = CONST_LINES[o + 2],
              x2 = CONST_LINES[o + 3], y2 = CONST_LINES[o + 4], z2 = CONST_LINES[o + 5];
        if(fwd[0]*x1 + fwd[1]*y1 + fwd[2]*z1 <= 0) continue;
        if(fwd[0]*x2 + fwd[1]*y2 + fwd[2]*z2 <= 0) continue;
        const sx1 = centerX + scale*(rRight[0]*x1 + rRight[1]*y1 + rRight[2]*z1),
              sy1 = centerY + scale*(rUp[0]*x1 + rUp[1]*y1 + rUp[2]*z1),
              sx2 = centerX + scale*(rRight[0]*x2 + rRight[1]*y2 + rRight[2]*z2),
              sy2 = centerY + scale*(rUp[0]*x2 + rUp[1]*y2 + rUp[2]*z2);
        d += "M" + sx1.toFixed(1) + "," + sy1.toFixed(1) + "L" + sx2.toFixed(1) + "," + sy2.toFixed(1);
    }
    constLinesPath.setAttribute("d", d);
}

function renderStarPaths(){
    let pickedD = "", pickedColor = null;
    bands.forEach(b => { b._mainD = ""; b._outerD = ""; b._hiddenD = "" });
    for(const star of visibleStars){
        const b = star.b;
        if(hoveredStar && hoveredStar.j === star.j){
            pickedD += starPath(star.sx, star.sy, b.size, b.inner, b.points, b.rotation);
            pickedColor = b.color;
        }else if(hiddenStars.has(star.j)){
            if(showHiddenMode) b._hiddenD += starPath(star.sx, star.sy, b.size, b.inner, b.points, b.rotation);
        }else if(b.border && b.borderWidth > 0){
            const shrink = b.size > 0 ? (b.size - b.borderWidth)/b.size : 1;
            b._outerD += starPath(star.sx, star.sy, b.size, b.inner, b.points, b.rotation);
            b._mainD += starPath(star.sx, star.sy, b.size*shrink, b.inner*shrink, b.points, b.rotation);
        }else{
            b._mainD += starPath(star.sx, star.sy, b.size, b.inner, b.points, b.rotation);
        }
    }
    bands.forEach(b => {
        if(b._mainPath) b._mainPath.setAttribute("d", b._mainD);
        if(b._outerPath) b._outerPath.setAttribute("d", b._outerD);
        if(b._hiddenPath) b._hiddenPath.setAttribute("d", b._hiddenD);
    });
    if(starsPickedPath){
        starsPickedPath.setAttribute("d", pickedD);
        if(pickedColor) starsPickedPath.setAttribute("fill", pickedColor);
        starsPickedPath.setAttribute("fill-opacity", hoveredStar && hoveredStar.hidden ? .7 : .3);
    }
}

function nearestStar(x, y){
    let best = null, bestDist = 1/0;
    for(const star of visibleStars){
        const dx = star.sx - x, dy = star.sy - y, dist = dx*dx + dy*dy;
        if(dist < bestDist){
            bestDist = dist;
            best = star;
        }
    }
    return best;
}

function updateHover(x, y){
    let next = null;
    if(pickMode && x != null){
        const star = nearestStar(x, y);
        if(star) next = {j: star.j, hidden: hiddenStars.has(star.j)}
    }
    if((!next && !hoveredStar) || (next && hoveredStar && next.j === hoveredStar.j && next.hidden === hoveredStar.hidden)) return;
    hoveredStar = next;
    renderStarPaths();
}

const barSvg = document.getElementById("bar"),
      BAR_LEFT = 46, BAR_WIDTH = 44, BAR_TOP = 14, BAR_HEIGHT = 410, BAR_BOTTOM = BAR_TOP + BAR_HEIGHT,
      SVG_NS = "http://www.w3.org/2000/svg";

function magToY(mag){
    return BAR_TOP + (mag - magMin)/(magMax - magMin)*BAR_HEIGHT;
}

function yToMag(y){
    return magMin + (y - BAR_TOP)/BAR_HEIGHT*(magMax - magMin);
}

function clampBarY(y){
    return Math.max(BAR_TOP, Math.min(BAR_BOTTOM, y));
}

let draggedBoundIdx = -1;

function makeSvgEl(tag, attrs){
    const el = document.createElementNS(SVG_NS, tag);
    for(const k in attrs){
        el.setAttribute(k, attrs[k]);
    }
    return el;
}

function rebuildBar(){
    while(barSvg.firstChild) barSvg.removeChild(barSvg.firstChild);
    for(let i = 0; i < bands.length; i += 1){
        const y0 = clampBarY(magToY(magBounds[i])), y1 = clampBarY(magToY(magBounds[i + 1]));
        if(y1 <= y0) continue;
        barSvg.appendChild(makeSvgEl("rect", {
            x: BAR_LEFT, y: y0, width: BAR_WIDTH, height: y1 - y0,
            fill: bands[i].render ? BAND_COLORS[i%BAND_COLORS.length] : "#dddddd"
        }));
    }
    barSvg.appendChild(makeSvgEl("g", {id: "barDots"}));
    for(let i = 1; i < magBounds.length - 1; i += 1){
        const y = magToY(magBounds[i]);
        if(y < BAR_TOP || y > BAR_BOTTOM) continue;
        barSvg.appendChild(makeSvgEl("line", {
            x1: BAR_LEFT - 2, x2: BAR_LEFT + BAR_WIDTH + 2, y1: y, y2: y, stroke: "#333333", "stroke-width": 3
        }));
        const handle = makeSvgEl("rect", {x: BAR_LEFT, y: y - 6, width: BAR_WIDTH, height: 12, fill: "transparent"});
        handle.style.cursor = "ns-resize";
        handle.addEventListener("pointerdown", e => {
            draggedBoundIdx = i;
            e.preventDefault();
            barSvg.setPointerCapture(e.pointerId);
        });
        barSvg.appendChild(handle);
    }
    const tickLo = Math.ceil(magMin), tickHi = Math.floor(magMax), tickStep = tickHi - tickLo > 7 ? 2 : 1;
    for(let m = tickLo; m <= tickHi; m += 1){
        const y = magToY(m);
        barSvg.appendChild(makeSvgEl("line", {
            x1: BAR_LEFT, x2: BAR_LEFT + BAR_WIDTH, y1: y, y2: y, stroke: "#000000", "stroke-width": 2, "stroke-dasharray": "4 3"
        }));
        if((m - tickLo)%tickStep === 0){
            const label = makeSvgEl("text", {x: BAR_LEFT + BAR_WIDTH + 6, y: y + 3, "font-size": 9, fill: "#999999"});
            label.textContent = m;
            barSvg.appendChild(label)
        }
    }
    rebuildBarDots();
}

function rebuildBarDots(){
    const dotsGroup = document.getElementById("barDots");
    if(!dotsGroup) return;
    while(dotsGroup.firstChild) dotsGroup.removeChild(dotsGroup.firstChild);
    const goldenStep = 2/(1 + Math.sqrt(5)), sortedMags = visibleMags.slice().sort((a, b) => a - b);
    for(let i = 0; i < sortedMags.length; i += 1){
        const y = magToY(sortedMags[i]);
        if(y < BAR_TOP || y > BAR_BOTTOM) continue;
        const jitter = goldenStep*i%1;
        dotsGroup.appendChild(makeSvgEl("circle", {cx: BAR_LEFT + 4 + jitter*(BAR_WIDTH - 8), cy: y, r: 1.6, fill: "#33333388"}));
    }
}

barSvg.addEventListener("pointermove", e => {
    if(draggedBoundIdx < 0) return;
    const rect = barSvg.getBoundingClientRect(), mag = yToMag(e.clientY - rect.top),
          lo = magBounds[draggedBoundIdx - 1] + .05, hi = magBounds[draggedBoundIdx + 1] - .05;
    magBounds[draggedBoundIdx] = Math.max(lo, Math.min(hi, mag));
    rebuildBar();
    repositionBandRows();
    render();
});
window.addEventListener("pointerup", () => { draggedBoundIdx = -1 });
barSvg.addEventListener("wheel", e => {
    e.preventDefault();
    const rect = barSvg.getBoundingClientRect(),
          pivot = yToMag(Math.max(BAR_TOP, Math.min(BAR_BOTTOM, e.clientY - rect.top))),
          factor = Math.exp(.0015*e.deltaY);
    let lo = pivot - (pivot - magMin)*factor, hi = pivot + (magMax - pivot)*factor;
    lo = Math.max(magMinBase, lo);
    hi = Math.min(magMaxBase, hi);
    if(hi - lo >= .2){
        magMin = lo;
        magMax = hi;
        rebuildBar();
        repositionBandRows();
    }
}, {passive: false});
barSvg.addEventListener("dblclick", () => {
    magMin = magMinBase;
    magMax = magMaxBase;
    rebuildBar();
    repositionBandRows();
});

const bandControlsEl = document.getElementById("bandcontrols");

function makeRowBreak(){
    const el = document.createElement("span");
    el.className = "rowbreak";
    return el;
}

function lighten(hex, amount){
    const n = parseInt(hex.slice(1), 16), g = n >> 8 & 255, b = 255 & n,
          mix = c => Math.round(c + (255 - c)*amount);
    return "#" + [mix(n >> 16 & 255), mix(g), mix(b)].map(c => c.toString(16).padStart(2, "0")).join("");
}

function subLabel(base, sub){
    const span = document.createElement("span");
    span.append(base);
    const subEl = document.createElement("sub");
    subEl.textContent = sub;
    span.appendChild(subEl);
    span.append(": ");
    return span;
}

function makeColorInput(value, title){
    const el = document.createElement("input");
    el.type = "text";
    el.className = "colorbox";
    el.value = value;
    el.title = title;
    return el;
}

function rebuildBandControls(){
    bandControlsEl.innerHTML = "";
    bandControlsEl.style.height = BAR_BOTTOM + BAR_TOP + "px";
    while(liveStarsGroup.firstChild) liveStarsGroup.removeChild(liveStarsGroup.firstChild);
    bands.forEach((band, idx) => {
        band._outerPath = document.createElementNS(SVG_NS, "path");
        band._outerPath.setAttribute("fill", band.borderColor || "#000000");
        liveStarsGroup.appendChild(band._outerPath);

        band._mainPath = document.createElementNS(SVG_NS, "path");
        band._mainPath.setAttribute("fill", band.color);
        liveStarsGroup.appendChild(band._mainPath);

        band._hiddenPath = document.createElementNS(SVG_NS, "path");
        band._hiddenPath.setAttribute("fill", band.color);
        band._hiddenPath.setAttribute("fill-opacity", "0.35");
        liveStarsGroup.appendChild(band._hiddenPath);

        const row = document.createElement("div");
        row.className = "bandrow";
        row.style.background = lighten(band.render ? BAND_COLORS[idx%BAND_COLORS.length] : "#dddddd", .6);

        if(!band.render){
            row.classList.add("hiddenrow");
            row.appendChild(document.createTextNode("Hidden"));
            bandControlsEl.appendChild(row);
            band._row = row;
            return;
        }

        let skipInput, innerInput;
        if(band.skip == null) band.skip = band.size > 0 ? skipOfRatio(band.points, band.inner/band.size) : 1;
        const syncInner = () => {
            band.inner = band.size*ratioOfSkip(band.points, band.skip);
            if(innerInput) innerInput.value = +band.inner.toFixed(1);
        };

        row.appendChild(document.createTextNode("{"));
        const pointsInput = document.createElement("input");
        pointsInput.type = "number";
        pointsInput.value = band.points;
        pointsInput.min = 3;
        pointsInput.max = 20;
        pointsInput.step = 1;
        pointsInput.title = "n in the Schläfli symbol {n/k}";
        pointsInput.addEventListener("input", () => {
            band.points = Math.max(3, parseInt(pointsInput.value) || 3);
            band.skip = wrapSkip(band.skip, band.points);
            if(skipInput) skipInput.value = +band.skip.toFixed(3);
            syncInner();
            render();
        });
        row.appendChild(pointsInput);
        row.appendChild(document.createTextNode("/"));

        skipInput = document.createElement("input");
        skipInput.type = "number";
        skipInput.value = +band.skip.toFixed(3);
        skipInput.min = 0;
        skipInput.step = 1;
        skipInput.title = "k in the Schläfli symbol {n/k}";
        skipInput.addEventListener("input", () => {
            band.skip = wrapSkip(parseFloat(skipInput.value) || 0, band.points);
            skipInput.value = +band.skip.toFixed(3);
            syncInner();
            render()
        });
        row.appendChild(skipInput);
        row.appendChild(document.createTextNode("}"));
        row.appendChild(makeRowBreak());

        const roLabel = subLabel("r", "o");
        roLabel.prepend("(");
        row.appendChild(roLabel);
        const sizeInput = document.createElement("input");
        sizeInput.type = "number";
        sizeInput.value = band.size;
        sizeInput.min = 0;
        sizeInput.step = 2;
        sizeInput.title = "Outer Radius";
        sizeInput.addEventListener("input", () => {
            band.size = parseFloat(sizeInput.value) || 0;
            syncInner();
            render()
        });
        row.appendChild(sizeInput);
        row.appendChild(document.createTextNode("px "));

        row.appendChild(subLabel("r", "i"));
        innerInput = document.createElement("input");
        innerInput.type = "number";
        innerInput.value = +band.inner.toFixed(1);
        innerInput.min = 0;
        innerInput.step = 1;
        innerInput.title = "Inner Radius";
        innerInput.addEventListener("input", () => {
            band.inner = parseFloat(innerInput.value) || 0;
            if(band.size > 0) band.skip = skipOfRatio(band.points, band.inner/band.size);
            if(skipInput) skipInput.value = +band.skip.toFixed(3);
            render();
        });
        row.appendChild(innerInput);
        row.appendChild(document.createTextNode("px)"));
        row.appendChild(makeRowBreak());

        const rotationInput = document.createElement("input");
        rotationInput.type = "number";
        rotationInput.value = +band.rotation.toFixed(1);
        rotationInput.step = 1;
        rotationInput.title = "Star Rotation";
        rotationInput.addEventListener("input", () => {
            band.rotation = parseFloat(rotationInput.value) || 0;
            render();
        });
        row.appendChild(rotationInput);
        row.appendChild(document.createTextNode("°"));

        const colorInput = makeColorInput(band.color, "Star Color");
        colorInput.addEventListener("input", () => {
            band.color = colorInput.value;
            band._mainPath.setAttribute("fill", band.color);
            band._hiddenPath.setAttribute("fill", band.color);
        });
        row.appendChild(colorInput);

        const borderChk = document.createElement("input");
        borderChk.type = "checkbox";
        borderChk.title = "Enable border";
        borderChk.checked = !!band.border;
        const borderLabel = document.createElement("label");
        borderLabel.className = "borderlabel";
        borderLabel.appendChild(borderChk);
        borderLabel.append("Border");
        row.appendChild(borderLabel);
        row.appendChild(makeRowBreak());

        const borderRow = document.createElement("span");
        borderRow.className = "borderrow";
        borderRow.style.display = band.border ? "" : "none";
        borderRow.appendChild(document.createTextNode("Border: "));
        const borderWidthInput = document.createElement("input");
        borderWidthInput.type = "number";
        borderWidthInput.value = band.borderWidth;
        borderWidthInput.min = 0;
        borderWidthInput.step = 1;
        borderWidthInput.title = "Border Width";
        borderWidthInput.addEventListener("input", () => {
            band.borderWidth = Math.max(0, parseFloat(borderWidthInput.value) || 0);
            render();
        });
        borderRow.appendChild(borderWidthInput);
        borderRow.appendChild(document.createTextNode("px "));

        const borderColorInput = makeColorInput(band.borderColor || "#000000", "Border Color");
        borderColorInput.addEventListener("input", () => {
            band.borderColor = borderColorInput.value;
            band._outerPath.setAttribute("fill", band.borderColor);
        });
        borderRow.appendChild(borderColorInput);
        row.appendChild(borderRow);

        borderChk.addEventListener("change", () => {
            band.border = borderChk.checked;
            borderRow.style.display = band.border ? "" : "none";
            repositionBandRows();
            render();
        });

        const splitBtn = document.createElement("button");
        splitBtn.textContent = "Split";
        splitBtn.title = "split this band into two";
        splitBtn.className = "split";
        splitBtn.addEventListener("click", () => splitBand(idx));
        row.appendChild(splitBtn);

        const deleteBtn = document.createElement("button");
        deleteBtn.textContent = "×";
        deleteBtn.title = "delete this section";
        deleteBtn.className = "del";
        deleteBtn.disabled = bands.length <= 2;
        deleteBtn.addEventListener("click", () => deleteBand(idx));
        row.appendChild(deleteBtn);

        bandControlsEl.appendChild(row);
        band._row = row;
    });
    repositionBandRows();
}

function deleteBand(idx){
    if(bands.length <= 2) return;
    if(idx === bands.length - 1) magBounds.splice(idx, 1);
    else magBounds.splice(idx + 1, 1);
    bands.splice(idx, 1);
    rebuildAll();
}

const ROW_GAP = 4;

function isotonicNondecreasing(values){
    const runs = [];
    for(let i = 0; i < values.length; i += 1){
        let run = {value: values[i], weight: 1, start: i, end: i};
        while(runs.length && runs[runs.length - 1].value > run.value){
            const prev = runs.pop(), weight = prev.weight + run.weight;
            run = {value: (prev.value*prev.weight + run.value*run.weight)/weight, weight, start: prev.start, end: run.end};
        }
        runs.push(run);
    }
    const out = new Array(values.length);
    for(const run of runs){
        for(let i = run.start; i <= run.end; i += 1){
            out[i] = run.value;
        }
    }
    return out;
}

function layoutRows(centers, heights, gap){
    const n = centers.length, offset = new Array(n).fill(0);
    for(let i = 1; i < n; i += 1){
        offset[i] = offset[i - 1] + (heights[i - 1] + heights[i])/2 + gap;
    }
    return isotonicNondecreasing(centers.map((c, i) => c - offset[i])).map((v, i) => v + offset[i]);
}

function repositionBandRows(){
    const rows = [];
    bands.forEach((band, i) => {
        if(!band._row) return;
        const y0 = clampBarY(magToY(magBounds[i])), y1 = clampBarY(magToY(magBounds[i + 1]));
        if(y1 - y0 < 1){
            band._row.style.display = "none";
        }else{
            band._row.style.display = "flex";
            rows.push({row: band._row, center: (y0 + y1)/2});
        }
    });
    if(!rows.length) return;
    const heights = rows.map(r => r.row.getBoundingClientRect().height),
          tops = layoutRows(rows.map(r => r.center), heights, ROW_GAP);
    rows.forEach((r, i) => { r.row.style.top = tops[i] + "px" });
}

function splitBand(idx){
    const wasLast = idx === bands.length - 1, mid = (magBounds[idx] + magBounds[idx + 1])/2;
    magBounds.splice(idx + 1, 0, mid);
    bands.splice(idx + 1, 0, {
        render: bands[idx].render, points: bands[idx].points, skip: bands[idx].skip,
        size: bands[idx].size, inner: bands[idx].inner, rotation: bands[idx].rotation,
        color: bands[idx].color, border: bands[idx].border,
        borderWidth: bands[idx].borderWidth, borderColor: bands[idx].borderColor
    });
    if(wasLast) bands[idx].render = true;
    rebuildAll();
}

function rebuildAll(){
    rebuildBar();
    rebuildBandControls();
    render();
}

let isPanning = false, lastPanX = 0, lastPanY = 0, panStartX = 0, panStartY = 0;
let isRotating = false, rotatePrevAngle = 0;

function flagCoords(e){
    const rect = flagSvg.getBoundingClientRect();
    return [(e.clientX - rect.left)/rect.width*flagW, (e.clientY - rect.top)/rect.height*flagH];
}

function pointerAngleDeg(e){
    const rect = flagSvg.getBoundingClientRect(),
          cx = rect.left + centerX/flagW*rect.width, cy = rect.top + centerY/flagH*rect.height;
    return Math.atan2(e.clientY - cy, e.clientX - cx)*180/Math.PI;
}

function wrapDeg(deg){
    return ((deg + 180)%360 + 360)%360 - 180;
}

function updateFlagCursor(){
    flagSvg.style.cursor = pickMode || moveCenterMode ? "crosshair" : isPanning || isRotating || isMovingCenter ? "grabbing" : "grab";
}

flagSvg.addEventListener("pointerdown", e => {
    if(e.button === 0){
        if(moveCenterMode){
            isMovingCenter = true;
            const [x, y] = flagCoords(e);
            centerX = x;
            centerY = y;
            render();
        }else if(e.shiftKey){
            isRotating = true;
            rotatePrevAngle = pointerAngleDeg(e);
        }else{
            isPanning = true;
            lastPanX = panStartX = e.clientX;
            lastPanY = panStartY = e.clientY;
        }
        flagSvg.classList.add("dragging");
        flagSvg.setPointerCapture(e.pointerId);
        updateFlagCursor();
    }
});
flagSvg.addEventListener("pointermove", e => {
    if(isMovingCenter){
        const [x, y] = flagCoords(e);
        centerX = x;
        centerY = y;
        render();
    }
    if(isRotating){
        const angle = pointerAngleDeg(e), delta = wrapDeg(angle - rotatePrevAngle), sign = mirrorView ? 1 : -1;
        setRoll((parseFloat(rollNumber.value) || 0) + sign*delta);
        rotatePrevAngle = angle;
        render();
    }
    if(isPanning){
        const rect = flagSvg.getBoundingClientRect(), pxPerScreenPx = flagW/rect.width,
              dx = (e.clientX - lastPanX)*pxPerScreenPx, dy = (e.clientY - lastPanY)*pxPerScreenPx;
        lastPanX = e.clientX;
        lastPanY = e.clientY;
        const scale = pxPerUnit*zoom, {rRight, rUp} = rolledAxes(),
              signX = mirrorView ? -1 : 1;
        rotateBasis(rUp, signX*-dx/scale);
        rotateBasis(rRight, dy/scale);
        reorthonormalizeBasis();
        render();
    }
    if(pickMode){
        const [x, y] = flagCoords(e);
        updateHover(x, y);
    }
});
flagSvg.addEventListener("pointerup", e => {
    isPanning = false;
    isRotating = false;
    isMovingCenter = false;
    flagSvg.classList.remove("dragging");
    updateFlagCursor();
    if(pickMode && Math.hypot(e.clientX - panStartX, e.clientY - panStartY) <= 5){
        const star = nearestStar(...flagCoords(e));
        if(star){
            if(hiddenStars.has(star.j)) hiddenStars.delete(star.j); else hiddenStars.add(star.j);
            console.log(hiddenStars);
            hoveredStar = {j: star.j, hidden: hiddenStars.has(star.j)};
            updateShowHiddenLabel();
            renderStarPaths();
        }
    }
});
flagSvg.addEventListener("pointercancel", () => {
    isPanning = false;
    isRotating = false;
    isMovingCenter = false;
    flagSvg.classList.remove("dragging");
    updateFlagCursor();
});
flagSvg.addEventListener("pointerleave", () => {
    if(pickMode) updateHover(null, null);
});
flagSvg.addEventListener("wheel", e => {
    e.preventDefault();
    let logZoom = Math.log10(zoom) - 9e-4*e.deltaY;
    logZoom = Math.max(-1.4, Math.min(1.4, logZoom));
    setZoom(Math.pow(10, logZoom));
    render();
}, {passive: false});

const zoomNumber = document.getElementById("zoomNum"), rollNumber = document.getElementById("rollNum"),
      raNumber = document.getElementById("raNum"), decNumber = document.getElementById("decNum"),
      reversedChk = document.getElementById("reversedChk"),
      wNumber = document.getElementById("wNum"), hNumber = document.getElementById("hNum");

function onZoomNumber(){
    const pct = Math.max(4, Math.min(2500, parseFloat(zoomNumber.value) || 100));
    zoom = pct/100;
    zoomNumber.value = Math.round(100*zoom);
    render();
}

function onRollNumber(){
    const deg = Math.max(-180, Math.min(180, parseFloat(rollNumber.value) || 0));
    rollRad = deg*Math.PI/180;
    rollNumber.value = Math.round(deg);
    render();
}

function setZoom(z){
    zoom = z;
    zoomNumber.value = Math.round(100*z);
}

function setRoll(deg){
    deg = ((deg + 180)%360 + 360)%360 - 180;
    rollRad = deg*Math.PI/180;
    rollNumber.value = Math.round(deg);
}

function setReflect(val){
    reflect = val;
}

function setMirrorView(val){
    mirrorView = val;
    reversedChk.checked = val;
}

function applyPointing(){
    rotmat = rotmatFromCentroid(vecFromRaDec(parseFloat(raNumber.value) || 0, parseFloat(decNumber.value) || 0));
    render();
}

function onRaNumber(){
    const deg = ((parseFloat(raNumber.value) || 0)%360 + 360)%360;
    raNumber.value = +deg.toFixed(1);
    applyPointing();
}

function onDecNumber(){
    const deg = Math.max(-90, Math.min(90, parseFloat(decNumber.value) || 0));
    decNumber.value = +deg.toFixed(1);
    applyPointing();
}

zoomNumber.addEventListener("input", onZoomNumber);
rollNumber.addEventListener("input", onRollNumber);
raNumber.addEventListener("input", onRaNumber);
decNumber.addEventListener("input", onDecNumber);
reversedChk.addEventListener("change", () => { mirrorView = reversedChk.checked; render(); });

wNumber.addEventListener("input", () => {
    const w = parseFloat(wNumber.value);
    if(w > 0) hNumber.value = Math.round(w*flagH/flagW);
});
hNumber.addEventListener("input", () => {
    const h = parseFloat(hNumber.value);
    if(h > 0) wNumber.value = Math.round(h*flagW/flagH);
});

const constellationSelect = document.getElementById("constellation");

function rotmatFromCentroid(centroid){
    const fwd = normalize(centroid);
    let right = cross([0, 0, 1], fwd);
    if(Math.hypot(right[0], right[1], right[2]) < 1e-6) right = cross([0, 1, 0], fwd);
    right = normalize(right);
    const up = cross(fwd, right);
    return [[-right[0], -right[1], -right[2]], [-up[0], -up[1], -up[2]], fwd.slice()];
}

function vecFromRaDec(raDeg, decDeg){
    const ra = raDeg*Math.PI/180, dec = decDeg*Math.PI/180;
    return [Math.cos(ra)*Math.cos(dec), Math.sin(ra)*Math.cos(dec), Math.sin(dec)];
}

function raDecFromVec(v){
    let ra = Math.atan2(v[1], v[0])*180/Math.PI;
    if(ra < 0) ra += 360;
    const dec = Math.asin(Math.max(-1, Math.min(1, v[2])))*180/Math.PI;
    return [ra, dec];
}

function syncPointingFields(fwd){
    const [ra, dec] = raDecFromVec(fwd);
    raNumber.value = +ra.toFixed(1);
    decNumber.value = +dec.toFixed(1);
}

function onConstellationChange(abbr){
    if(abbr === flag.default_constellation){
        rotmat = flag.rotmat.map(v => v.slice());
        setReflect(flag.reflect);
        setMirrorView(flag.reflect);
        setRoll(flag.roll0);
    }else{
        rotmat = rotmatForConstellation(abbr);
        setReflect(false);
        setMirrorView(false);
        setRoll(0);
    }
    render();
}

function rotmatForConstellation(abbr){
    return rotmatFromCentroid(CONSTELLATIONS.find(c => c.abbr === abbr).c);
}

const flagSelect = document.getElementById("flagsel");

const COUNTRY_KEYS = ["AU", "BR", "NZ", "PNG", "WS"];
const SUBDIVISION_COUNTRY = {
    SANTACRUZ: "Argentina", TDF: "Argentina",
    ACT: "Australia", NSW: "Australia", NT: "Australia", VIC: "Australia", CC: "Australia", CX: "Australia",
    MAG: "Chile",
    NIUE: "New Zealand", TOK: "New Zealand",
    SONSOROL: "Palau",
    ENB: "Papua New Guinea", NIRE: "Papua New Guinea", SIMBU: "Papua New Guinea",
    WNB: "Papua New Guinea", WPNG: "Papua New Guinea",
    AK: "United States",
};
const COUNTRY_ABBR = {
    Argentina: "AR", Australia: "AU", Chile: "CL", "New Zealand": "NZ",
    Palau: "PW", "Papua New Guinea": "PG", "United States": "US",
};

function stripRedundantCountry(name, country){
    const suffix = ` (${country})`;
    return name.endsWith(suffix) ? name.slice(0, -suffix.length) : name;
}

function buildFlagOptions(){
    if(FLAGS.Blank){
        const opt = document.createElement("option");
        opt.value = "Blank";
        opt.textContent = "Blank";
        flagSelect.appendChild(opt);
    }
    const byName = (a, b) => FLAGS[a].name.localeCompare(FLAGS[b].name);
    const countries = FLAG_ORDER.filter(k => COUNTRY_KEYS.includes(k)).sort(byName);
    const subdivisions = FLAG_ORDER.filter(k => k in SUBDIVISION_COUNTRY).sort((a, b) => {
        const c = SUBDIVISION_COUNTRY[a].localeCompare(SUBDIVISION_COUNTRY[b]);
        return c !== 0 ? c : byName(a, b);
    });
    const used = new Set([...countries, ...subdivisions, "Blank"]);
    const others = FLAG_ORDER.filter(key => !used.has(key)).sort(byName);

    const addGroup = (label, keys, labelFor) => {
        if(!keys.length) return;
        const group = document.createElement("optgroup");
        group.label = label;
        keys.forEach(key => {
            const opt = document.createElement("option");
            opt.value = key;
            opt.textContent = labelFor ? labelFor(key) : FLAGS[key].name;
            group.appendChild(opt);
        });
        flagSelect.appendChild(group);
    };
    addGroup("Countries", countries);
    addGroup("Subdivisions", subdivisions, key => {
        const country = SUBDIVISION_COUNTRY[key];
        return `${COUNTRY_ABBR[country] || country}/${stripRedundantCountry(FLAGS[key].name, country)}`;
    });
    addGroup("Other", others);

    return countries[0] || subdivisions[0] || others[0];
}

function selectFlag(key){
    fetch("svg/flag_clips/"+key+".svg")
        .then(response => response.text())
        .then(text => {
            flag = FLAGS[key];
            flagW = flag.W;
            flagH = flag.H;
            flagSvg.setAttribute("viewBox", `${flag.minX || 0} ${flag.minY || 0} ${flagW} ${flagH}`);
            flagSvg.style.height = 760*flagH/flagW + "px";

            let defs = "", clipAttr = "";

            const doc = new DOMParser().parseFromString(text, "image/svg+xml");
            const shape = doc.documentElement.innerHTML;
            console.log("svg/flag_clips/"+key+".svg");
            console.log(shape);

            defs = `<mask id="liveMask">${shape}</mask>`;
            //defs = `<clipPath id="liveClip">${shape}</clipPath>`;
            clipAttr = ' mask="url(#liveMask)"';
            flagSvg.innerHTML = flag.base + defs + `<g${clipAttr}><path id="constLines"></path><g id="liveStarsGroup"></g><path id="liveStarsHi"></path></g>` +
                `<g id="centerCrosshair" style="display:none"><path id="crosshairOutline" fill="none" stroke="#000000" stroke-width="8" vector-effect="non-scaling-stroke"></path><path id="crosshairFg" fill="none" stroke="#ee0000" stroke-width="4" vector-effect="non-scaling-stroke"></path></g>`;

            constLinesPath = document.getElementById("constLines");
            liveStarsGroup = document.getElementById("liveStarsGroup");
            starsPickedPath = document.getElementById("liveStarsHi");
            centerCrosshairG = document.getElementById("centerCrosshair");
            crosshairOutline = document.getElementById("crosshairOutline");
            crosshairFg = document.getElementById("crosshairFg");
            constLinesPath.setAttribute("fill", "none");
            constLinesPath.setAttribute("stroke-width", "2");
            constLinesPath.setAttribute("vector-effect", "non-scaling-stroke");
            centerCrosshairG.style.display = moveCenterMode ? "" : "none";

            const isBlank = key === "Blank";
            bgColorLabel.style.display = isBlank ? "" : "none";
            bgColorInput.style.display = isBlank ? "" : "none";
            if(isBlank){
                bgColorInput.value = "#000000";
                const bgEl = document.getElementById("blankBg");
                if(bgEl) bgEl.setAttribute("fill", bgColorInput.value);
            }
            updateLineColor();

            if(flag.rotmat === null){
                rotmat = rotmatForConstellation(flag.default_constellation);
            }else{
                rotmat = flag.rotmat.map(v => v.slice());
            }
            setReflect(flag.reflect);
            setMirrorView(flag.reflect);
            pxPerUnit = flag.scale;
            centerX = flag.cx;
            centerY = flag.cy;
            clipRegion = flag.region.slice();
            magBounds = flag.bounds.slice();
            bands = flag.bands.map(b => ({
                render: b.render, points: b.points, skip: b.skip, size: b.size, inner: b.inner, rotation: b.rotation,
                color: b.color || "#ffffff", border: !!b.border,
                borderWidth: b.borderWidth || 0, borderColor: b.borderColor || "#000000"
            }));
            magMinBase = flag.magMin;
            magMaxBase = flag.magMax;
            magMin = magMinBase;
            magMax = magMaxBase;
            hiddenStars = new Set(flag.hidden || []);
            updateShowHiddenLabel();
            hoveredStar = null;

            setZoom(1);
            setRoll(flag.roll0);
            wNumber.value = 1920;
            hNumber.value = Math.round(1920*flagH/flagW);
            constellationSelect.value = flag.default_constellation;
            updateFlagCursor();
            rebuildBar();
            rebuildBandControls();
            render();
            updateOriginalView();
        });
    // if(flag.clip === "circle"){
    //     clipAttr = ' clip-path="url(#Bstars)"'
    // }else{
    //     const response = await fetch("svg/flag_clips/"+key+".svg");
    //     const text = await response.text();
    //     const doc = new DOMParser().parseFromString(text, "image/svg+xml");

    //     let shape = doc.documentElement.innerHTML;

    //     defs = `<clipPath id="liveClip">${shape}</clipPath>`;
    //     clipAttr = ' clip-path="url(#liveClip)"'
    // }
}

const SVG_ORIGINAL_DIR = "svg/original";

const showOriginalChk = document.getElementById("showOriginalChk"),
      origFlagImg = document.getElementById("origFlagImg"),
      liveControls = document.getElementById("liveControls");

function updateOriginalView(){
    const showOriginal = showOriginalChk.checked;
    flagSvg.style.display = showOriginal ? "none" : "";
    origFlagImg.style.display = showOriginal ? "" : "none";
    if(showOriginal) origFlagImg.src = encodeURI(SVG_ORIGINAL_DIR + "/" + flag.originalFile);
    liveControls.classList.toggle("disabled", showOriginal);
    for(const el of liveControls.querySelectorAll("input, select, button")) el.disabled = showOriginal;
}
showOriginalChk.addEventListener("change", updateOriginalView);

const bgColorInput = document.getElementById("bgColorInput"),
      bgColorLabel = document.getElementById("bgColorLabel");

function relativeLuminance(hex){
    const n = parseInt(hex.slice(1), 16), r = n >> 16 & 255, g = n >> 8 & 255, b = 255 & n;
    return (0.299*r + 0.587*g + 0.114*b)/255;
}

function updateLineColor(){
    if(!constLinesPath) return;
    const bg = bgColorInput.style.display !== "none" ? bgColorInput.value : (flag && flag.bg),
          closeToWhite = bg && relativeLuminance(bg) >= .9;
    constLinesPath.setAttribute("stroke", closeToWhite ? "#000000" : "#ffffff");
}

bgColorInput.addEventListener("input", () => {
    const bgEl = document.getElementById("blankBg");
    if(bgEl) bgEl.setAttribute("fill", bgColorInput.value);
    updateLineColor();
});

function exportSvgClone(outW, outH){
    const pngW = Math.round(outW || Math.max(flagW, 1200)), pngH = Math.round(outH || pngW*flagH/flagW),
          clone = flagSvg.cloneNode(true);
    clone.setAttribute("width", pngW);
    clone.setAttribute("height", pngH);
    clone.style.background = "";
    clone.style.width = "";
    clone.style.height = "";
    clone.style.display = "";
    const cloneCrosshair = clone.querySelector("#centerCrosshair");
    if(cloneCrosshair) cloneCrosshair.remove();
    return {svgText: new XMLSerializer().serializeToString(clone), pngW, pngH};
}

function triggerDownload(blob, filename){
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    URL.revokeObjectURL(a.href);
}

async function fetchOriginalSvgText(){
    const response = await fetch(encodeURI(SVG_ORIGINAL_DIR + "/" + flag.originalFile));
    return await response.text();
}

document.getElementById("downloadSvgBtn").addEventListener("click", async () => {
    const filenameBase = flag.name.replace(/[^\w]+/g, "_");
    if(showOriginalChk.checked){
        const svgText = await fetchOriginalSvgText();
        triggerDownload(new Blob([svgText], {type: "image/svg+xml"}), filenameBase + "_original.svg");
    }else{
        const {svgText} = exportSvgClone();
        triggerDownload(new Blob([svgText], {type: "image/svg+xml"}), filenameBase + "_stars.svg");
    }
});

document.getElementById("downloadPngBtn").addEventListener("click", async () => {
    const filenameBase = flag.name.replace(/[^\w]+/g, "_"),
          outW = parseFloat(wNumber.value) || Math.max(flagW, 1200),
          outH = parseFloat(hNumber.value) || outW*flagH/flagW,
          pngW = Math.round(outW), pngH = Math.round(outH);
    try {
        let svgText;
        if(showOriginalChk.checked){
            const doc = new DOMParser().parseFromString(await fetchOriginalSvgText(), "image/svg+xml");
            doc.documentElement.setAttribute("width", pngW);
            doc.documentElement.setAttribute("height", pngH);
            svgText = new XMLSerializer().serializeToString(doc);
        }else{
            ({svgText} = exportSvgClone(outW, outH));
        }
        const blob = await flag_to_png(svgText, pngW, pngH);
        triggerDownload(blob, filenameBase + (showOriginalChk.checked ? "_original.png" : "_stars.png"));
    } catch(err) {
        console.error(err);
        alert("Sorry — this browser couldn't rasterize the flag to PNG.");
    }
});

const hideBtn = document.getElementById("hideBtn"), showBtn = document.getElementById("showBtn"),
      showHiddenLabel = document.getElementById("showHiddenLabel"),
      moveCenterBtn = document.getElementById("moveCenterBtn");

function updateShowHiddenLabel(){
    showHiddenLabel.textContent = `Show ${hiddenStars.size} Hidden`;
}

function setPickMode(val){
    pickMode = val;
    hideBtn.classList.toggle("active", pickMode);
    if(pickMode) setMoveCenterMode(false);
    updateFlagCursor();
    if(!pickMode) updateHover(null, null);
}
hideBtn.addEventListener("click", () => setPickMode(!pickMode));

function setMoveCenterMode(val){
    moveCenterMode = val;
    moveCenterBtn.classList.toggle("active", moveCenterMode);
    if(centerCrosshairG) centerCrosshairG.style.display = moveCenterMode ? "" : "none";
    if(moveCenterMode) setPickMode(false);
    updateFlagCursor();
}
moveCenterBtn.addEventListener("click", () => setMoveCenterMode(!moveCenterMode));
window.addEventListener("keydown", e => {
    if(e.key === "Escape" && pickMode) setPickMode(false);
    if(e.key === "Escape" && moveCenterMode) setMoveCenterMode(false);
});
showBtn.addEventListener("change", () => {
    showHiddenMode = showBtn.checked;
    renderStarPaths();
});

const linesBtn = document.getElementById("linesBtn");
linesBtn.addEventListener("change", () => {
    showConstLines = linesBtn.checked;
    render();
});

const onlyLinesBtn = document.getElementById("onlyLinesBtn");
onlyLinesBtn.addEventListener("change", () => {
    onlyLinesMode = onlyLinesBtn.checked;
    render();
});

async function init(){
    const data = await fetch("flag_star_data.json").then(r => r.json());
    FLAGS = data.flags;
    FLAG_ORDER = data.FLAG_ORDER;
    CONSTELLATIONS = data.CONSTELLATIONS;
    BAND_COLORS = data.BAND_COLORS;
    STAR_DATA = new Float32Array(data.STAR_DATA);
    CONST_LINES = new Float32Array(data.CONST_LINES);
    CONST_LINE_RANGES = data.CONST_LINE_RANGES || {};

    CONSTELLATIONS.forEach(c => {
        const opt = document.createElement("option");
        opt.value = c.abbr;
        opt.textContent = c.name;
        constellationSelect.appendChild(opt);
    });
    constellationSelect.addEventListener("change", () => onConstellationChange(constellationSelect.value));

    const firstKey = buildFlagOptions();
    flagSelect.addEventListener("change", () => selectFlag(flagSelect.value));
    document.getElementById("reset").addEventListener("click", () => selectFlag(flagSelect.value));

    flagSelect.value = firstKey;
    selectFlag(firstKey);
}

init();
