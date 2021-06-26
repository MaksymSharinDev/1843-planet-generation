/*
 * From https://www.redblobgames.com/x/1843-planet-generation/
 * Copyright 2018 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 *
 * Adapting mapgen4 code for a sphere. Quick & dirty, for procjam2018
 */
'use strict';

const SimplexNoise = require('simplex-noise');
const FlatQueue = require('flatqueue');
const {vec3, mat4} = require('gl-matrix');
const {makeRandInt, makeRandFloat} = require('@redblobgames/prng');
const SphereMesh = require('./sphere-mesh');
const maps = require('./maps');
const renderEngine = require('./render');

const WATER_LEVEL = 0;

/* UI parameters */
let N = 10000;
let P = 20;
let jitter = 0.75;
let rotation = -1;
let tilt = 0.75;
let procession = 1.1;
let drawMode = 'centroid';
let draw_plateVectors = false;
let draw_plateBoundaries = false;
let draw_landBoundaries = false;
let draw_elevationLines = false;
let draw_equator = false;
let draw_extraLat = true;
let draw_primeMeridian = false;
let draw_windVectors = false;
let draw_normalVectors = false;

let renderLayer = "surface";

let SEA_LEVEL = 0.5;
let SEED = 123; // 41 is pretty good too

let mapNeedsRedraw = true;
let globeNeedsRedraw = true;
let mapProjectionType = maps.equirectangular;
let draw_windVectors_map = true;

window.setN = newN => { N = newN; generateMesh(); };
window.setP = newP => { P = newP; generateMap(); };
window.setSeed = newSeed => { SEED = newSeed; generateMap(); };
window.setSeaLevel = newLevel => { SEA_LEVEL = newLevel; generateMap(); };
window.setJitter = newJitter => { jitter = newJitter; generateMesh(); };
window.setRotation = newRotation => { rotation = newRotation; draw(); };
window.setTilt     = newTilt     => { tilt = newTilt; draw(); };
window.setProcession = newProcession => { procession = newProcession; draw(); };
window.setDrawMode = newMode => { drawMode = newMode; draw(); };
window.setDrawPlateVectors = flag => { draw_plateVectors = flag; draw(); };
window.setDrawPlateBoundaries = flag => { draw_plateBoundaries = flag; draw(); };
window.setDrawLandBoundaries = flag => { draw_landBoundaries = flag; draw(); }
window.setDrawElevationLines = flag => { draw_elevationLines = flag; draw(); }
window.setDrawEquator = flag => { draw_equator = flag; draw(); };
window.setDrawExtraLat = flag => { draw_extraLat = flag; draw(); };
window.setDrawPrimeMeridian = flag => { draw_primeMeridian = flag; draw(); };
window.setDrawExtraLat = flag => { draw_extraLat = flag; draw(); };
window.setDrawWindVectors = flag => { draw_windVectors = flag; draw(); };
window.setDrawNormalVectors = flag => { draw_normalVectors = flag; draw(); };

window.setRender = renderWhat => { renderLayer = renderWhat; draw(); };
window.advanceWeather = () => { advanceWeather(mesh, map); draw(); }

let advanceWeatherIntervalId;
window.startWeather = () => {
    advanceWeatherIntervalId = window.setInterval(() => {
        advanceWeather(mesh, map); 
        draw();
    }, 10);
};
window.pauseWeather = () => {
    clearInterval(advanceWeatherIntervalId); 
};

window.setMapProjection = projectionName => {
    if (projectionName === "mercator") mapProjectionType = maps.mercator;
    else if (projectionName === "equirectangular") mapProjectionType = maps.equirectangular;
    else if (projectionName === "sinusoidal") mapProjectionType = maps.sinusoidal;
    else if (projectionName === "naturalearth") mapProjectionType = maps.natural_earth;

    mapNeedsRedraw = true;
    draw();
};
window.setDrawWindVectors_map = flag => {draw_windVectors_map = flag; mapNeedsRedraw = true; draw(); }

/**********************************************************************
 * Geometry
 */

let _randomNoise = new SimplexNoise(makeRandFloat(SEED));
const persistence = 2/3;
const amplitudes = Array.from({length: 5}, (_, octave) => Math.pow(persistence, octave));

function fbm_noise(nx, ny, nz) {
    let sum = 0, sumOfAmplitudes = 0;
    for (let octave = 0; octave < amplitudes.length; octave++) {
        let frequency = 1 << octave;
        sum += amplitudes[octave] * _randomNoise.noise3D(nx * frequency, ny * frequency, nz * frequency);
        sumOfAmplitudes += amplitudes[octave];
    }
    return sum / sumOfAmplitudes;
}

/* Calculate the centroid and push it onto an array */
function pushCentroidOfTriangle(out, ax, ay, az, bx, by, bz, cx, cy, cz) {
    // TODO: renormalize to radius 1
    out.push((ax+bx+cx)/3, (ay+by+cy)/3, (az+bz+cz)/3);
}


function generateTriangleCenters(mesh, {r_xyz}) {
    let {numTriangles} = mesh;
    let t_xyz = [];
    for (let t = 0; t < numTriangles; t++) {
        let a = mesh.s_begin_r(3*t),
            b = mesh.s_begin_r(3*t+1),
            c = mesh.s_begin_r(3*t+2);
        pushCentroidOfTriangle(t_xyz,
                 r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
                 r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
                 r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]);
    }
    return t_xyz;
}

function generateVoronoiGeometry(mesh, {r_xyz, t_xyz}, r_color_fn) {
    const {numSides} = mesh;
    let xyz = [], tm = [];

    for (let s = 0; s < numSides; s++) {
        let inner_t = mesh.s_inner_t(s),
            outer_t = mesh.s_outer_t(s),
            begin_r = mesh.s_begin_r(s);
        let rgb = r_color_fn(begin_r);
        xyz.push(t_xyz[3*inner_t], t_xyz[3*inner_t+1], t_xyz[3*inner_t+2],
                      t_xyz[3*outer_t], t_xyz[3*outer_t+1], t_xyz[3*outer_t+2],
                      r_xyz[3*begin_r], r_xyz[3*begin_r+1], r_xyz[3*begin_r+2]);
        tm.push(rgb, rgb, rgb);
    }
    return {xyz, tm};
}

class QuadGeometry {
    constructor () {
        /* xyz = position in 3-space;
           tm = temperature, moisture
           I = indices for indexed drawing mode */
    }

    setMesh({numSides, numRegions, numTriangles}) {
        this.I = new Int32Array(3 * numSides);
        this.xyz = new Float32Array(3 * (numRegions + numTriangles));
        this.tm = new Float32Array(2 * (numRegions + numTriangles));
    }

    setMap(mesh, {r_xyz, t_xyz, r_color_fn, s_flow, r_elevation, t_elevation, r_moisture, t_moisture}) {
        const V = 0.95;
        const {numSides, numRegions, numTriangles} = mesh;
        const {xyz, tm, I} = this;

        xyz.set(r_xyz);
        xyz.set(t_xyz, r_xyz.length);
        // TODO: multiply all the r, t points by the elevation, taking V into account

        let p = 0;
        for (let r = 0; r < numRegions; r++) {
            tm[p++] = r_elevation[r];
            tm[p++] = r_moisture[r];
        }
        for (let t = 0; t < numTriangles; t++) {
            tm[p++] = t_elevation[t];
            tm[p++] = t_moisture[t];
        }

        let i = 0, count_valley = 0, count_ridge = 0;
        let {_halfedges, _triangles} = mesh;
        for (let s = 0; s < numSides; s++) {
            let opposite_s = mesh.s_opposite_s(s),
                r1 = mesh.s_begin_r(s),
                r2 = mesh.s_begin_r(opposite_s),
                t1 = mesh.s_inner_t(s),
                t2 = mesh.s_inner_t(opposite_s);
            
            // Each quadrilateral is turned into two triangles, so each
            // half-edge gets turned into one. There are two ways to fold
            // a quadrilateral. This is usually a nuisance but in this
            // case it's a feature. See the explanation here
            // https://www.redblobgames.com/x/1725-procedural-elevation/#rendering
            let coast = r_elevation[r1] < 0.0 || r_elevation[r2] < 0.0;
            if (coast || s_flow[s] > 0 || s_flow[opposite_s] > 0) {
                // It's a coastal or river edge, forming a valley
                I[i++] = r1; I[i++] = numRegions+t2; I[i++] = numRegions+t1;
                count_valley++;
            } else {
                // It's a ridge
                I[i++] = r1; I[i++] = r2; I[i++] = numRegions+t1;
                count_ridge++;
            }
        }

        console.log('ridge=', count_ridge, ', valley=', count_valley);
    }
}

/**********************************************************************
 * Plates
 */

function pickRandomRegions(mesh, N, randInt) {
    let {numRegions} = mesh;
    let chosen_r = new Set();
    while (chosen_r.size < N && chosen_r.size < numRegions) {
        chosen_r.add(randInt(numRegions));
    }
    return chosen_r;
}


function generatePlates(mesh, r_xyz) {
    let r_plate = new Int32Array(mesh.numRegions);
    r_plate.fill(-1);
    let plate_r = pickRandomRegions(mesh, Math.min(P, N), makeRandInt(SEED));
    let queue = Array.from(plate_r);
    for (let r of queue) { r_plate[r] = r; }
    let out_r = [];
    const randInt = makeRandInt(SEED);

    /* In Breadth First Search (BFS) the queue will be all elements in
       queue[queue_out ... queue.length-1]. Pushing onto the queue
       adds an element to the end, increasing queue.length. Popping
       from the queue removes an element from the beginning by
       increasing queue_out.

       To add variety, use a random search instead of a breadth first
       search. The frontier of elements to be expanded is still
       queue[queue_out ... queue.length-1], but pick a random element
       to pop instead of the earliest one. Do this by swapping
       queue[pos] and queue[queue_out].
    */
    
    for (let queue_out = 0; queue_out < queue.length; queue_out++) {
        let pos = queue_out + randInt(queue.length - queue_out);
        let current_r = queue[pos];
        queue[pos] = queue[queue_out];
        mesh.r_circulate_r(out_r, current_r);
        for (let neighbor_r of out_r) {
            if (r_plate[neighbor_r] === -1) {
                r_plate[neighbor_r] = r_plate[current_r];
                queue.push(neighbor_r);
            }
        }
    }

    // Assign a random movement vector for each plate
    let plate_vec = [];
    for (let center_r of plate_r) {
        let neighbor_r = mesh.r_circulate_r([], center_r)[0];
        let p0 = r_xyz.slice(3 * center_r, 3 * center_r + 3),
            p1 = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
        plate_vec[center_r] = vec3.normalize([], vec3.subtract([], p1, p0));
    }

    return {plate_r, r_plate, plate_vec};
}


/* Distance from any point in seeds_r to all other points, but 
 * don't go past any point in stop_r */
function assignDistanceField(mesh, seeds_r, stop_r) {
    const randInt = makeRandInt(SEED);
    let {numRegions} = mesh;
    let r_distance = new Float32Array(numRegions);
    r_distance.fill(Infinity);
    
    let queue = [];
    for (let r of seeds_r) {
        queue.push(r);
        r_distance[r] = 0;
    }

    /* Random search adapted from breadth first search */
    let out_r = [];
    for (let queue_out = 0; queue_out < queue.length; queue_out++) {
        let pos = queue_out + randInt(queue.length - queue_out);
        let current_r = queue[pos];
        queue[pos] = queue[queue_out];
        mesh.r_circulate_r(out_r, current_r);
        for (let neighbor_r of out_r) {
            if (r_distance[neighbor_r] === Infinity && !stop_r.has(neighbor_r)) {
                r_distance[neighbor_r] = r_distance[current_r] + 1;
                queue.push(neighbor_r);
            }
        }
    }
    return r_distance;
    // TODO: possible enhancement: keep track of which seed is closest
    // to this point, so that we can assign variable mountain/ocean
    // elevation to each seed instead of them always being +1/-1
}


/* Calculate the collision measure, which is the amount
 * that any neighbor's plate vector is pushing against 
 * the current plate vector. */
const COLLISION_THRESHOLD = 0.75;
function findCollisions(mesh, r_xyz, plate_is_ocean, r_plate, plate_vec) {
    const deltaTime = 1e-2; // simulate movement
    let {numRegions} = mesh;
    let mountain_r = new Set(),
        coastline_r = new Set(),
        ocean_r = new Set();
    let r_out = [];
    /* For each region, I want to know how much it's being compressed
       into an adjacent region. The "compression" is the change in
       distance as the two regions move. I'm looking for the adjacent
       region from a different plate that pushes most into this one*/
    for (let current_r = 0; current_r < numRegions; current_r++) {
        let bestCompression = Infinity, best_r = -1;
        mesh.r_circulate_r(r_out, current_r);
        for (let neighbor_r of r_out) {
            if (r_plate[current_r] !== r_plate[neighbor_r]) {
                /* sometimes I regret storing xyz in a compact array... */
                let current_pos = r_xyz.slice(3 * current_r, 3 * current_r + 3),
                    neighbor_pos = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
                /* simulate movement for deltaTime seconds */
                let distanceBefore = vec3.distance(current_pos, neighbor_pos),
                    distanceAfter = vec3.distance(vec3.add([], current_pos, vec3.scale([], plate_vec[r_plate[current_r]], deltaTime)),
                                                  vec3.add([], neighbor_pos, vec3.scale([], plate_vec[r_plate[neighbor_r]], deltaTime)));
                /* how much closer did these regions get to each other? */
                let compression = distanceBefore - distanceAfter;
                /* keep track of the adjacent region that gets closest */
                if (compression < bestCompression) {
                    best_r = neighbor_r;
                    bestCompression = compression;
                }
            }
        }
        if (best_r !== -1) {
            /* at this point, bestCompression tells us how much closer
               we are getting to the region that's pushing into us the most */
            let collided = bestCompression > COLLISION_THRESHOLD * deltaTime;
            if (plate_is_ocean.has(current_r) && plate_is_ocean.has(best_r)) {
                (collided? coastline_r : ocean_r).add(current_r);
                // ocean_r.add(current_r);
            } else if (!plate_is_ocean.has(current_r) && !plate_is_ocean.has(best_r)) {
                if (collided) mountain_r.add(current_r);
            } else {
                (collided? mountain_r : coastline_r).add(current_r);
            }
        }

        // TODO: divergent plate boundaries - trenches and canyon seas
    }

    // TODO: trench_r for divergent oceanic boundaries
    return {mountain_r, coastline_r, ocean_r};
}

// TODO: make hotspot trails (hawaii and such) by picking a starting point and following
// the opposite direction of plate movement. make the end a wider blob but also submerged
// TODO: make trails of shallow ocean starting from continents and extending in opposite direction of plate motion
// TODO: make ocean-ocean convergence make trench on the side of the larger plate and coastline
// on the side of the smaller
// TODO: make ocean-continent convergence create trench and mountains
// TODO: when doing these checks, make sure to check region type, not whole plate type (most plates are combination ocean/continent)
// REF: Earth's REAL Lost Continents by Atlas Pro on Youtube
function assignRegionElevation(mesh, {r_xyz, plate_is_ocean, r_plate, plate_vec, /* out */ r_elevation}) {
    const epsilon = 1e-3;
    let {numRegions} = mesh;

    let {mountain_r, coastline_r, ocean_r} = findCollisions(
        mesh, r_xyz, plate_is_ocean, r_plate, plate_vec);

    for (let r = 0; r < numRegions; r++) {
        if (r_plate[r] === r) {
            (plate_is_ocean.has(r)? ocean_r : coastline_r).add(r);
        }
    }

    let stop_r = new Set();
    for (let r of mountain_r) { stop_r.add(r); }
    for (let r of coastline_r) { stop_r.add(r); }
    for (let r of ocean_r) { stop_r.add(r); }

    console.log('seeds mountain/coastline/ocean:', mountain_r.size, coastline_r.size, ocean_r.size, 'plate_is_ocean', plate_is_ocean.size,'/', P);
    let r_distance_a = assignDistanceField(mesh, mountain_r, ocean_r);
    let r_distance_b = assignDistanceField(mesh, ocean_r, coastline_r);
    let r_distance_c = assignDistanceField(mesh, coastline_r, stop_r);

    for (let r = 0; r < numRegions; r++) {
        let a = r_distance_a[r] + epsilon,
            b = r_distance_b[r] + epsilon,
            c = r_distance_c[r] + epsilon;
        if (a === Infinity && b === Infinity) {
            r_elevation[r] = 0.1;
        } else {
            r_elevation[r] = (1/a - 1/b) / (1/a + 1/b + 1/c);
        }
        r_elevation[r] += 0.1 * fbm_noise(r_xyz[3*r], r_xyz[3*r+1], r_xyz[3*r+2]);
        
        r_elevation[r] -= SEA_LEVEL;
    }
}

function reigonIsWater(r_elevation, r) {
    return r_elevation[r] < WATER_LEVEL;
}



/**********************************************************************
 * Rivers - from mapgen4
 */

function assignTriangleValues(mesh, {r_elevation, r_moisture, /* out */ t_elevation, t_moisture}) {
    const {numTriangles} = mesh;
    for (let t = 0; t < numTriangles; t++) {
        let s0 = 3*t;
        let r1 = mesh.s_begin_r(s0),
            r2 = mesh.s_begin_r(s0+1),
            r3 = mesh.s_begin_r(s0+2);
        t_elevation[t] = 1/3 * (r_elevation[r1] + r_elevation[r2] + r_elevation[r3]);
        t_moisture[t] = 1/3 * (r_moisture[r1] + r_moisture[r2] + r_moisture[r3]);
    }
}


let _queue = new FlatQueue();
function assignDownflow(mesh, {t_elevation, /* out */ t_downflow_s, /* out */ order_t}) {
    /* Use a priority queue, starting with the ocean triangles and
     * moving upwards using elevation as the priority, to visit all
     * the land triangles */
    let {numTriangles} = mesh,
        queue_in = 0;
    t_downflow_s.fill(-999);
    /* Part 1: ocean triangles get downslope assigned to the lowest neighbor */
    for (let t = 0; t < numTriangles; t++) {
        if (t_elevation[t] < 0) {
            let best_s = -1, best_e = t_elevation[t];
            for (let j = 0; j < 3; j++) {
                let s = 3 * t + j,
                    e = t_elevation[mesh.s_outer_t(s)];
                if (e < best_e) {
                    best_e = e;
                    best_s = s;
                }
            }
            order_t[queue_in++] = t;
            t_downflow_s[t] = best_s;
            _queue.push(t, t_elevation[t]);
        }
    }
    /* Part 2: land triangles get visited in elevation priority */
    for (let queue_out = 0; queue_out < numTriangles; queue_out++) {
        let current_t = _queue.pop();
        for (let j = 0; j < 3; j++) {
            let s = 3 * current_t + j;
            let neighbor_t = mesh.s_outer_t(s); // uphill from current_t
            if (t_downflow_s[neighbor_t] === -999 && t_elevation[neighbor_t] >= 0.0) {
                t_downflow_s[neighbor_t] = mesh.s_opposite_s(s);
                order_t[queue_in++] = neighbor_t;
                _queue.push(neighbor_t, t_elevation[neighbor_t]);
            }
        }
    }
}


function assignFlow(mesh, {order_t, t_elevation, t_moisture, t_downflow_s, /* out */ t_flow, /* out */ s_flow}) {
    let {numTriangles, _halfedges} = mesh;
    s_flow.fill(0);
    for (let t = 0; t < numTriangles; t++) {
        if (t_elevation[t] >= 0.0) {
            t_flow[t] = 0.5 * t_moisture[t] * t_moisture[t];
        } else {
            t_flow[t] = 0;
        }
    }
    for (let i = order_t.length-1; i >= 0; i--) {
        let tributary_t = order_t[i];
        let flow_s = t_downflow_s[tributary_t];
        let trunk_t = (_halfedges[flow_s] / 3) | 0;
        if (flow_s >= 0) {
            t_flow[trunk_t] += t_flow[tributary_t];
            s_flow[flow_s] += t_flow[tributary_t]; // TODO: isn't s_flow[flow_s] === t_flow[?]
            if (t_elevation[trunk_t] > t_elevation[tributary_t]) {
                t_elevation[trunk_t] = t_elevation[tributary_t];
            }
        }
    }
}


/**********************************************************************
 * Weaather
 */

function vectorSubtract(a, b) {
    let [ax, ay, az] = a,
        [bx, by, bz] = b;
    
    return [ax - bx, ay - by, az - bz];
}

function dot (a, b) {
    let [ax, ay, az] = a,
        [bx, by, bz] = b;
    
    return ax*bx + ay*by + az*bz;
}

function magnitude(a) {
    return Math.sqrt(dot(a, a));
}

function getNextNeighbor(mesh, current_r, dir, {r_xyz}) {
    let bestAngle = Math.Infinity;
    let bestNeighbor = -1;

    let current_xyz = r_xyz.slice(3 * current_r, 3 * current_r + 3);

    let r_out = [];
    mesh.r_circulate_r(r_out, current_r);
    for (let neighbor_r of r_out) {
        let neighbor_xyz = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
        let neighbor_dir = vectorSubtract(neighbor_xyz, current_xyz);

        let angleProportionate = dot(dir, neighbor_dir) / magnitude(neighbor_dir);

        if (bestAngle === Math.Infinity || angleProportionate < bestAngle) {
            bestAngle = angleProportionate;
            bestNeighbor = neighbor_r;
        }
    }

    return bestNeighbor;
}

function getNeighbor(mesh, current_r, dir, {r_xyz}) {
    let dist = magnitude(dir);
    let current_xyz = r_xyz.slice(3 * current_r, 3 * current_r + 3);
    let n_count = 0;

    while (n_count++ < 20) {
        if (dist <= 0) return current_r;

        let next_r = getNextNeighbor(mesh, current_r, dir, {r_xyz});
        let next_xyz = r_xyz.slice(3 * next_r, 3 * next_r + 3);
        dist -= magnitude(vectorSubtract(next_xyz, current_xyz));

        current_r = next_r;
        current_xyz = next_xyz;
    }

    console.error("Vector was too long, or something else went wrong.");
    console.log({dir, mag: magnitude(dir)});

    return current_r;
}

function statsAnalysis_vectorMagnitude(vectors) {
    let data = [];
    for (let r = 0; r < vectors.length; r++) {
        data.push(magnitude(vectors[r]));
    }
    statsAnalysis(data);
}
function statsAnalysis_neighborsDistance(mesh, {r_xyz}) {
    let {numRegions} = mesh;
    let distances = [];
    for (let r = 0; r < numRegions; r++) {
        let current_xyz = r_xyz.slice(3 * r, 3 * r + 3);
        let r_out = [];
        mesh.r_circulate_r(r_out, r);
        for (let neighbor_r of r_out) {
            let neighbor_xyz = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
            let neighbor_dir = vectorSubtract(neighbor_xyz, current_xyz);

            distances.push(magnitude(neighbor_dir));
        }
    }
    statsAnalysis(data);
}
function statsAnalysis(data, do_histogram=false) {
    let distances = data;

    let average = 0;
    let max = -1;
    let min = 99999999;
    for (let i = 0; i < distances.length; i++) {
        average += distances[i];
        max = Math.max(max, distances[i]);
        min = Math.min(min, distances[i]);
    }
    average /= distances.length;
    
    let stddev = 0;
    for (let i = 0; i < distances.length; i++) {
        let tempx = distances[i] - average;
        stddev += tempx*tempx;
    }
    stddev = Math.sqrt(stddev / distances.length);

    let bins, stepsize;
    if(do_histogram) {
        const binCount = 50;
        stepsize = (max-min) / binCount;
        bins = new Array(binCount);
        bins.fill(0);
        for (let i = 0; i < distances.length; i++) {
            let d = distances[i];
            d -= min;
            let bin = Math.floor(d / stepsize);
            bins[bin]++;
        }
    }

    console.log({min, max, average, stddev, bins: bins+"", bins_stepsize:stepsize});
}

function xyzToLatLon(xyz, planetRadius = 1) {
    let [x, y, z] = xyz;
    let lat_deg = (180/Math.PI) * Math.acos(Math.abs(z) / planetRadius), 
        lon_deg = (180/Math.PI) * Math.atan2(y, x);
    let abs_lat_deg = Math.abs(lat_deg);

    return [lat_deg, lon_deg, abs_lat_deg];
}

function xyzFromLatLon([lat, lon], planetRadius = 1) {
    let x = planetRadius * Math.cos(lon) * Math.sin(lat),
        y = planetRadius * Math.sin(lon) * Math.sin(lat),
        z = planetRadius * Math.cos(lat);
    return [x, y, z];
}

function sigmoid(x) {
    return 1 / (1 + Math.exp(-x));
}

// note: the distance between neighboring tiles follows this distribution:
// { min: 0.00008017334191617021, max: 0.08736101006422599, average: 0.04041248913501134, stddev: 0.013437966505337578 }
// currrently, wind speeds follow this:
// { min: 0, max: 0.09183723630842007, average: 0.01448420366929869, stddev: 0.00795346069581971 }
function assignRegionWindVectors(mesh, {r_xyz, r_elevation, r_temperature, /* out */ r_wind}) {
    const planetRadius = 1;
    let {numRegions} = mesh;

    let temp_temp_winddirs_0 = [];
    let temp_temp_winddirs_1 = [];

    // TODO: try normalizing wind_dir just before converting to xyz (ie before calculating d_lat and d_lon)
    // The smaller its magnitude, the more tangent the final vector will be to the sphere
    // and thus the more accurate the final vector will be

    // Note that the wind vectors will always look weird unless viewed straight on
    // because they're lines tangent to the sphere, NOT lines laying flat on it like a shadow

    // TODO: store r_latlon in the map data
    // TODO: store r_wind_latlon ~~alongside~~ instead of r_wind_xyz
    // TODO: make versions of getNextNeighbor and friends that accept a latlon vector instead of an xyz
    //      NOTE: these functions will need to account for 357 being closer to 11 than 45 is to 11.

    const TEMPERATURE_INFLUENCE_FACTOR = 4;

    for (let r = 0; r < numRegions; r++) {
        let wind_speed = 0.02;
        // let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        // let lat_deg = (180/Math.PI) * Math.acos(Math.abs(z) / planetRadius), 
        //     lon_deg = (180/Math.PI) * Math.atan2(y, x);
        // let abs_lat_deg = Math.abs(lat_deg);
        let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        let [lat_deg, lon_deg, abs_lat_deg] = xyzToLatLon([x, y, z]);
        if (z < 0) z = -z;

        let wind_dir = [0, 0];

        // prevailing winds
        if (0 < abs_lat_deg && abs_lat_deg < 30) {
            let trigterm = (Math.PI * (lat_deg-0 )) / (2 * 30);
            wind_dir = [-Math.sin(trigterm), Math.cos(trigterm)];
        } else
        if (30 < abs_lat_deg && abs_lat_deg < 60) {
            let trigterm = (Math.PI * (lat_deg-30)) / (2 * 30);
            wind_dir = [Math.cos(trigterm), -Math.sin(trigterm)];
        } else 
        if (60 < abs_lat_deg && abs_lat_deg < 90) {
            let trigterm = (Math.PI * (lat_deg-60)) / (2 * 30);
            wind_dir = [-Math.sin(trigterm), Math.cos(trigterm)];
        }

        if (z < 0) {
            // southern hemisphere
            wind_dir[1] = -wind_dir[1];
        }

        // slowdown/speedup according to elevation change
        let blowsPast_r = getNextNeighbor(mesh, r, xyzFromLatLon(wind_dir), map);
        let elevation_change = Math.max(r_elevation[r], WATER_LEVEL) - Math.max(r_elevation[blowsPast_r], WATER_LEVEL);
        wind_speed += 0.1*elevation_change;//*= 2*(sigmoid(elevation_change)-0.5); //*= (sigmoid(elevation_change)+0.5);
        if(isNaN(wind_speed)) {
            console.log(sigmoid(elevation_change)+0.5);
            crash
        }

        // Wind blows from cold to warm
        // TODO: this causes NaN temperatures and NaN wind
        if (!isNaN(r_temperature[r])) {
            let r_out = [];
            mesh.r_circulate_r(r_out, r);
            for (let neighbor_r of r_out) {
                if (isNaN(r_temperature[neighbor_r])) continue;
                let [nx, ny, nz] = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
                let [n_lat_deg, n_lon_deg, n_abs_lat_deg] = xyzToLatLon([nx, ny, nz]);

                
                let d_lat = n_lat_deg - lat_deg,
                    d_lon = n_lon_deg - lon_deg;
                let mag = magnitude([d_lat, d_lon, 0]);
                if (mag === 0) continue;

                let temperatureDifference = r_temperature[neighbor_r] - r_temperature[r];
                let influence = TEMPERATURE_INFLUENCE_FACTOR * temperatureDifference;
                let influenceFactor = influence / mag;

                wind_dir[0] += d_lat * influenceFactor;
                wind_dir[1] += d_lon * influenceFactor;

            }

            // add speed according to how much temperature difference there is and how in-line the wind is with that temp diff
            wind_speed += (0.01 / TEMPERATURE_INFLUENCE_FACTOR) * Math.sqrt(wind_dir[0]*wind_dir[0] + wind_dir[1]*wind_dir[1]);
        }

        // theta is the around, phi is the up and down
        // theta is longitude, phi is lattitude
        let [dtheta, dphi] = wind_dir;
        temp_temp_winddirs_0.push(wind_speed*wind_dir[0]);
        temp_temp_winddirs_1.push(wind_speed*wind_dir[1]);
        let wind_blowsTo_lat_deg = lat_deg + dphi,
            wind_blowsTo_lon_deg = lon_deg + dtheta;
        let wind_blowsTo_lat_rad = (Math.PI/180) * wind_blowsTo_lat_deg,
            wind_blowsTo_lon_rad = (Math.PI/180) * wind_blowsTo_lon_deg;

        
        let wind_blowsTo_x = planetRadius * Math.cos(wind_blowsTo_lon_rad) * Math.sin(wind_blowsTo_lat_rad),
            wind_blowsTo_y = planetRadius * Math.sin(wind_blowsTo_lon_rad) * Math.sin(wind_blowsTo_lat_rad),
            wind_blowsTo_z = planetRadius * Math.cos(wind_blowsTo_lat_rad);

        let wind_x = wind_blowsTo_x - x,
            wind_y = wind_blowsTo_y - y,
            wind_z = wind_blowsTo_z - z;
        
        let mag = magnitude([wind_x, wind_y, wind_z]);
        mag = mag === 0 ? 1 : mag;
        let speedFactor = wind_speed / mag;
        r_wind[r] = [speedFactor*wind_x, speedFactor*wind_y, speedFactor*wind_z];

        // if (z > 0) r_wind[r] = [x,y,z];

        // // proof that the wind vectors are indeed tangent to the sphere
        // // even though they draw funny:
        // let dot = wind_x * x + wind_y * y + wind_z * z;
        // let angle = (180/Math.PI) * Math.acos(dot/100);
        // if (angle > 92 || angle < 88)
        // {
        //     console.log(angle);
        //     crash;
        // }

        // if (dot([wind_x, wind_y, wind_z], [x, y, z]) != -0.0602189418300005) {
        //     console.log("something went wrong");
        //     console.log(dot([wind_x, wind_y, wind_z], [x, y, z]))
            
        //     let wind_blowsTo_cartesian = [wind_blowsTo_x, wind_blowsTo_y, wind_blowsTo_z];
        //     let wind_blowsTo_spherical = [wind_blowsTo_lat_deg, wind_blowsTo_lon_deg];
        //     let wind = [wind_x, wind_y, wind_z];
        //     let r_loc= [x, y, z];
        //     console.log({lat_deg, lon_deg, dtheta, dphi, r_loc, wind_blowsTo_cartesian, wind_blowsTo_spherical, wind});
            
            
        //     crash;
        // }
    }

    // statsAnalysis(temp_temp_winddirs_0);
    // statsAnalysis(temp_temp_winddirs_1);
    // statsAnalysis_vectorMagnitude(r_wind, numRegions);
    

    // statsAnalysis_vectorMagnitude(r_wind, numRegions);
}

function assignRegionTemperature(mesh, {r_xyz, r_elevation, /* out */ r_temperature}) {
    const planetRadius = 1;
    let {numRegions} = mesh;

    for (let r = 0; r < numRegions; r++) {
        let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        let lat_deg = (180/Math.PI) * Math.acos(Math.abs(z) / planetRadius), 
            lon_deg = (180/Math.PI) * Math.atan2(y, x);

        r_temperature[r] = (1-r_elevation[r]) * lat_deg / 90;
    }
}

function assignRegionHumidity(mesh, {r_elevation, r_temperature, /* out */ r_humidity}) {
    let {numRegions} = mesh;

    for (let r = 0; r < numRegions; r++) {
        r_humidity[r] = reigonIsWater(r_elevation, r) ? 0.5 * r_temperature[r] : 0;
    }
}

function assignRegionClouds(mesh, {/* out */ r_clouds}) {
    let {numRegions} = mesh;
    for (let r = 0; r < numRegions; r++) {
        r_clouds[r] = 0;
    }
}

function average(list) {
    return sum(list) / list.length;
}

function sum(list) {
    let retval = 0;
    for (let i = 0; i < list.length; i++) {
        retval += list[i];
    }
    return retval;
}

// TODO: go back to the way the reference does it
function reassignRegionTemperature(mesh, {r_xyz, r_elevation, r_wind, r_clouds, /* in/out */ r_temperature}) {
    const planetRadius = 1;
    let {numRegions} = mesh;
    let r_newTemperature = new Array(numRegions);
    
    let r_blownFrom = new Array(numRegions);

    const PERSISTENCE = 0.85;
    
    for (let r = 0; r < numRegions; r++) {
        let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        let lat_deg = (180/Math.PI) * Math.acos(Math.abs(z) / planetRadius), 
            lon_deg = (180/Math.PI) * Math.atan2(y, x);

        // base temperature
        let baseTemperature = (1-r_elevation[r]) * lat_deg / 90;
        baseTemperature = PERSISTENCE*r_temperature[r] + (1-PERSISTENCE)*baseTemperature;

        // account for cloud cover
        if (isNaN(r_clouds[r])) r_clouds[r] = 0; // I have no idea how r_clouds[r] can ever be NaN, so this line is my bandaid
        baseTemperature -= 2 * (1-Math.max(1, 2*r_clouds[r]));

        // average with old temperature
        // r_newTemperature[r] = 0.25 * baseTemperature + 0.75 * r_temperature[r];
        r_newTemperature[r] = baseTemperature;

        // wind
        let blownTo_r;
        try {
            blownTo_r = getNeighbor(mesh, r, r_wind[r], {r_xyz});    
        } catch (error) {
            console.error(error);
            blownTo_r = r;
        }

        if (!r_blownFrom[blownTo_r]) r_blownFrom[blownTo_r] = [];
        r_blownFrom[blownTo_r].push(r); 
    }

    // diffusion
    for (let r = 0; r < numRegions; r++) {
        let r_out = [];
        mesh.r_circulate_r(r_out, r);
        let neighborAverage = 0;
        for (let neighbor_r of r_out) {
            neighborAverage += r_temperature[neighbor_r];
        }
        neighborAverage /= r_out.length;
        if (r_out.length == 0)  console.error("no neighbors");

        r_newTemperature[r] = (3/4) * r_newTemperature[r] + (1/4) * neighborAverage;
    }
    
    // account for rain
    for(let r = 0; r < numRegions; r++) {
        // rain occurs above cloud value 0.5
        r_newTemperature[r] -= 2 * Math.max(0, r_clouds[r]-0.5);
    }

    // account for wind
    for(let r = 0; r < numRegions; r++) {
        // temperature blown in from neighbors
        let blownTemperatures = r_blownFrom[r] 
            ? r_blownFrom[r].map(neighbor_r => r_newTemperature[neighbor_r])
            : [r_newTemperature[r]];
        if (blownTemperatures.length === 0) blownTemperatures = [r_newTemperature[r]];

        let tempBlown = average(blownTemperatures);
        r_newTemperature[r] = (r_newTemperature[r] + tempBlown) / 2;

        // cooldown from wind blowing
        r_newTemperature[r] -= 2*(magnitude(r_wind[r]) - 0.01);
    }

    // finalize
    // TODO: try putting humidity and temperature through a sigmoid funciton
    for (let r = 0; r < numRegions; r++) {
        r_temperature[r] = r_newTemperature[r]; //Math.max(0, 0.5*(0.5*r_baseTemperature[r] + 0.5*r_newTemperature[r]));
    }
}

function reassignRegionHumidity(mesh, {r_xyz, r_elevation, r_wind, r_clouds, r_temperature, /* in/out */ r_humidity}) {
    const planetRadius = 1;
    let {numRegions} = mesh;
    let r_newHumidity = new Array(numRegions);
    let blownHumidity = new Array(numRegions);
    blownHumidity.fill([]);
    
    const PERSISTENCE = 0.8;

    for (let r = 0; r < numRegions; r++) {
        // let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        // let lat_deg = (180/Math.PI) * Math.acos(Math.abs(z) / planetRadius), 
        //     lon_deg = (180/Math.PI) * Math.atan2(y, x);

        // base humidity
        let baseHumidity = reigonIsWater(r_elevation, r) ? 0.5 * r_temperature[r] : 0;

        // // average with old humidity
        // r_newHumidity[r] = 1.25 * baseHumidity + 0.5 * r_humidity[r];

        r_newHumidity[r] = PERSISTENCE*r_humidity[r] + (1-PERSISTENCE)*baseHumidity;

        // wind
        let blownTo_r;
        try {
            blownTo_r = getNeighbor(mesh, r, r_wind[r], {r_xyz});    
        } catch (error) {
            console.error(error);
            blownTo_r = r;
        }
        
        if (!blownHumidity[blownTo_r]) blownHumidity[blownTo_r] = [];
        blownHumidity[blownTo_r].push(r_humidity[r]); 
    }

    // account for rain
    for(let r = 0; r < numRegions; r++) {
        // rain occurs above cloud value 0.5
        r_newHumidity[r] -= 0.01*2*Math.max(0, r_clouds[r]-0.5);
    }

    // account for wind
    // TODO: do the humidity blown, humidity kept thing from clouds
    for(let r = 0; r < numRegions; r++) {
        let humidityBlown = average(blownHumidity[r] ? blownHumidity[r] : [r_newHumidity[r]]);
        // if (isNaN(blownHumidity[r])) blownHumidity[r] = r_humidity[r]; // some tiles don't get wind blowing onto them

        r_newHumidity[r] = (r_newHumidity[r] + humidityBlown) / 2;
        // r_temperature[r] += 0.5*(magnitude(r_wind[r]) - 5);
    }
    
    // diffusion
    for (let r = 0; r < numRegions; r++) {
        let r_out = [];
        mesh.r_circulate_r(r_out, r);
        let neighborAverage = 0;
        for (let neighbor_r of r_out) {
            neighborAverage += r_newHumidity[neighbor_r];
        }
        neighborAverage /= r_out.length;

        r_newHumidity[r] = (3/4) * r_newHumidity[r] + (1/4) * neighborAverage;
    }

    // finalize
    // TODO: try putting humidity and temperature through a sigmoid funciton
    for (let r = 0; r < numRegions; r++) {
        r_humidity[r] = Math.max(0, r_newHumidity[r]);
    }
}

function reassignRegionClouds(mesh, {r_xyz, r_wind, r_temperature, r_humidity, /* out */ r_clouds}) {
    const planetRadius = 1;
    let {numRegions} = mesh;
    let r_newClouds = new Array(numRegions);
    r_newClouds.fill(0);
    let blownClouds = new Array(numRegions);
    
    let cloudChallenge = 2;
    let cloudFactorCalc = x => (x + 1) * Math.pow(x, cloudChallenge) / 2;

    for (let r = 0; r < numRegions; r++) {
        // new clouds
        // note: this is not the method used by the reference; I wanted to have a continuous cloud system instead of a discrete one
        //r_newClouds[r] += (1-r_temperature[r]) * r_humidity[r];
        r_newClouds[r] = cloudFactorCalc(1-r_temperature[r]) * cloudFactorCalc(r_humidity[r]);

        // keep existing clouds, but discount some dissapation rate
        if (isNaN(r_clouds[r])) r_clouds[r] = 0; // I have no idea how r_clouds[r] can ever be NaN, so this line is my bandaid
        r_newClouds[r] += 0.8*r_clouds[r];

        // if it's raining, decrease the cloud level
        // note: it's raining when clouds are above 0.5
        r_newClouds[r] -= 2*Math.max(0, r_clouds[r] - 0.5);

        // wind
        let blownTo_r;
        try {
            blownTo_r = getNeighbor(mesh, r, r_wind[r], {r_xyz});    
        } catch (error) {
            console.error(error);
            blownTo_r = r;
        }
        
        let windspeed = magnitude(r_wind[r]);
        let cloudsBlown =r_newClouds[r]; //Math.max(0, Math.min(1, windspeed / 0.04)) * r_newClouds[r];
        let cloudsKept = 0;//(1-Math.max(0, Math.min(1, windspeed / 0.04))) * r_newClouds[r];
        r_newClouds[r] = cloudsKept;
        
        if (isNaN(r_newClouds[r])) {
            console.log({cloudsBlown, cloudsKept, windspeed, temperature: r_temperature[r], humidity: r_humidity[r]})
            crash
        }

        if (!blownClouds[blownTo_r]) blownClouds[blownTo_r] = [];
        blownClouds[blownTo_r].push(cloudsBlown); 
    }
    
    for (let r = 0; r < numRegions; r++) {
        let cloudsBlown = sum(blownClouds[r] ? blownClouds[r] : [r_newClouds[r]]);
        // if (isNaN(blownClouds[r])) blownClouds[r] = r_clouds[r];
        r_newClouds[r] += cloudsBlown;
    }

    
    // diffusion and finalize
    for (let r = 0; r < numRegions; r++) {
        let r_out = [];
        mesh.r_circulate_r(r_out, r);
        let neighborAverage = 0;
        for (let neighbor_r of r_out) {
            neighborAverage += r_newClouds[neighbor_r];
        }
        neighborAverage /= r_out.length;

        r_clouds[r] = (1/4) * r_newClouds[r] + (3/4) * neighborAverage;
        r_clouds[r] = Math.min(1, Math.max(0, r_clouds[r]));
    }
    
    // // finalize
    // for (let r = 0; r < numRegions; r++) {
    //     r_clouds[r] = r_newClouds[r];
    // }

    // statsAnalysis(r_clouds)
}

function advanceWeather(mesh, map) {
    assignRegionWindVectors(mesh, map);
    reassignRegionTemperature(mesh, map);
    reassignRegionHumidity(mesh, map);
    reassignRegionClouds(mesh, map);

    globeNeedsRedraw = true;

    // console.log("Stats for temperature, humidity, and clouds");
    // statsAnalysis(map.r_temperature);
    // statsAnalysis(map.r_humidity);
    // statsAnalysis(map.r_clouds);
}



/**********************************************************************
 * Main
 */

// ugh globals, sorry
var mesh, map = {};
var quadGeometry = new QuadGeometry();

function generateMesh() {
    let result = SphereMesh.makeSphere(N, jitter, makeRandFloat(SEED));
    mesh = result.mesh;
    quadGeometry.setMesh(mesh);
    
    map.r_xyz = result.r_xyz;
    map.t_xyz = generateTriangleCenters(mesh, map);
    map.r_elevation = new Float32Array(mesh.numRegions);
    map.t_elevation = new Float32Array(mesh.numTriangles);
    map.r_moisture = new Float32Array(mesh.numRegions);
    map.t_moisture = new Float32Array(mesh.numTriangles);
    map.t_downflow_s = new Int32Array(mesh.numTriangles);
    map.order_t = new Int32Array(mesh.numTriangles);
    map.t_flow = new Float32Array(mesh.numTriangles);
    map.s_flow = new Float32Array(mesh.numSides);
    map.r_wind = new Array(mesh.numRegions);
    map.r_temperature = new Array(mesh.numRegions);
    map.r_humidity = new Array(mesh.numRegions);
    map.r_clouds = new Array(mesh.numRegions);
    
    generateMap();

    // TODO: eventually, simulate a game of Slay to draw country borders
}

function generateMap() {
    let result = generatePlates(mesh, map.r_xyz);
    map.plate_r = result.plate_r;
    map.r_plate = result.r_plate;
    map.plate_vec = result.plate_vec;
    map.plate_is_ocean = new Set();
    for (let r of map.plate_r) {
        if (makeRandInt(r)(10) < 5) {
            map.plate_is_ocean.add(r);
            // TODO: either make tiny plates non-ocean, or make sure tiny plates don't create seeds for rivers
        }
    }
    assignRegionElevation(mesh, map);
    // TODO: assign region moisture in a better way!
    // I'll be adding in a version of this method: https://nickmcd.me/2018/07/10/procedural-weather-patterns/
    // Combined with a simple prevailing winds model and major currents model
    for (let r = 0; r < mesh.numRegions; r++) {
        map.r_moisture[r] = (map.r_plate[r] % 10) / 10.0;
    }
    assignTriangleValues(mesh, map);
    assignDownflow(mesh, map);
    assignFlow(mesh, map);

    assignRegionClouds(mesh, map);
    assignRegionWindVectors(mesh, map);
    assignRegionTemperature(mesh, map);
    assignRegionHumidity(mesh, map);

    quadGeometry.setMap(mesh, map);
    draw();
}

let _draw_pending = false;
function _draw() {
    renderEngine.draw({
        all: {
            mesh,
            map,
            N,

            layer: renderLayer,
            generateVoronoiGeometry: generateVoronoiGeometry,

            draw_axis: true,

            draw_plateVectors,
            draw_plateBoundaries,
            draw_windVectors,
            draw_normalVectors,
            draw_equator,
            draw_extraLat,
            draw_primeMeridian,
            draw_landBoundaries,
            draw_elevationLines,
        },

        globe: {
            forceRedraw: globeNeedsRedraw,
            surfaceMode: drawMode,
            procession,
            tilt,
            rotation: -rotation,

            cloud_height: 0.1,
        },

        map: {
            forceRedraw: mapNeedsRedraw,
            draw_windVectors: draw_windVectors_map,
            projection: mapProjectionType,

            draw_axis: false,
        },
    });
    mapNeedsRedraw = false;
    globeNeedsRedraw = false;

    _draw_pending = false;
}

function draw() {
    if (!_draw_pending) {
        _draw_pending = true;
        requestAnimationFrame(_draw);
    }
}

generateMesh();
