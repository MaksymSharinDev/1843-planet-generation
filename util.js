
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

// the two functions below took a TON of fiddling, but I think I finally have them right
function xyzFromLatLon_rad(latRad, lonRad) {
    // https://math.stackexchange.com/a/989911
    return [Math.cos(latRad) * Math.cos(lonRad),
            Math.sin(latRad),
            Math.cos(latRad) * Math.sin(lonRad)];
}

function xyzToLatLon_rad(x, y, z) {
    // formulas from https://gis.stackexchange.com/a/120685
    let r = magnitude([x,y,z]);
    if (r === 0) return undefined;

    let atan2 = (z, x) => {
        if (x === 0) {
            if (z > 0) return  Math.PI; // should be pi/2 ?
            if (z < 0) return -Math.PI; // should be pi/2 ?
            console.log("x and z was zero")
            return undefined;
        }

        if (z === 0) {
            if (x < 0) return  Math.PI; // should be pi/2 ?
            if (x > 0) return -Math.PI; // should be pi/2 ?
            console.log("z and x was zero")
            return undefined;
        } 

        if (x >= 0) {
            return Math.atan(z / x);// + Math.PI/2;
        }
        // something that should be returning PI is returning PI/2 instead

        if (z >= 0) {
            // if (Math.pow(Math.atan(z/x), 2) < 2 - Math.pow(Math.asin(y / r),2)) return Math.PI; // this is the case where things go sideways
            // there's some set of x,z such that Math.atan(z/x) ~= sqrt(2 - Math.pow(Math.asin(y / r),2))) when it should ~= pi/2
            // find that set

            // let expval = Math.sqrt(2 - Math.pow(Math.asin(y / r),2));
            // let explat = Math.abs(Math.atan(z/x));
            // let epsilon = 0.1;
            // if (expval - epsilon <= explat && explat <= expval + epsilon) return 1.5*Math.PI;

            // do 1.5*Math.PI to view the map, leave it as Math.PI for accurate wind vectors
            return Math.atan(z/x) + Math.PI; //2*Math.PI; // this *2 arguably makes it worse, but lets lines render for some reason
        } else {
            return Math.atan(z/x) - Math.PI;
        }
    }

    return [
        Math.asin(y / r),
        atan2(z, x)//Math.atan2(z, x) // issues when x is neg and z is pos
    ];
}

const DEG2RAD = Math.PI / 180;
const RAD2DEG = 180 / Math.PI;

function xyzFromLatLon_deg(latDeg, lonDeg) {
    return xyzFromLatLon_rad(latDeg * DEG2RAD, lonDeg * DEG2RAD);
}

function xyzToLatLon_deg(x, y, z) {
    let [lat_rad, lon_rad] = xyzToLatLon_rad(x, y, z);
    return [lat_rad * RAD2DEG, lon_rad * RAD2DEG];
}

function normalize(xyz) {
    let [x, y, z] = xyz;
    let mag = magnitude(xyz);

    return [x/mag, y/mag, z/mag];
}

function test_xyzlatlon() {
    // let latlon = [
    //     [0, 0],
    //     [10, 0],
    //     [20, 0],
    //     [-10, 0],
    //     [-20, 0],

    //     [90, 0],
    //     [-90, 0],
    //     [180, 0], // any lattitude greater than 180 or less than -180 is invalid, but these functions handle these cases very gracefully
    //     [-180, 0], // and interpret them generously to mean [0, 180]. beautiful :)
        
    //     [0, 10],
    //     [0, -10],
    //     [0, 360],
    //     [0, 370],
    // ];

    // let testU = ([lat, lon]) => {
    //     let [x, y, z] = xyzFromLatLon_deg(lat, lon);

    //     console.log([lat, lon] + " -> " + [x, y, z] + " -> " + xyzToLatLon_deg(x, y, z));
    // }

    // latlon.forEach(e => testU(e));

    console.log(xyzToLatLon_deg(-1, 0, 0));
    console.log(xyzToLatLon_deg(1,  0, 0));
    console.log(xyzToLatLon_deg(0, 0, -1));
    console.log(xyzToLatLon_deg(0, 0,  1));
    
    
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

function clamp(min, max, val) {
    return Math.min(max, Math.max(min, val));
}

/*
    Takes:
        a [delta_lattitude, delta_longitude] vector (in degrees),
        the location on the sphere [lattitude, longitude] (in degrees) where the vector's origin is,
        and the xyz location corresponding to the latlon location

    Returns the vector in cartesian coordinates
*/
function deg_latlonVectorAt_to_xyzVector(latlonVec_deg, latlonAt_deg, [x, y, z]) {
    let [lat_deg, lon_deg] = latlonAt_deg;
    let [dtheta, dphi] = latlonVec_deg;
    
    let pointsTo_lat_deg = lat_deg + dphi,
        pointsTo_lon_deg = lon_deg + dtheta;
    let [pointsTo_x, pointsTo_y, pointsTo_z] = xyzFromLatLon_deg(pointsTo_lat_deg, pointsTo_lon_deg);

    return [
        pointsTo_x - x,
        pointsTo_y - y,
        pointsTo_z - z
    ];
}

function setMagnitude(xyzVector, targetMagnitude) {
    let [x, y, z] = xyzVector;
    let mag = magnitude(xyzVector);
    let scale = targetMagnitude / mag;
    return [scale*x, scale*y, scale*z];
}

module.exports = { 
    clamp, 

    statsAnalysis, 
    statsAnalysis_neighborsDistance, 
    statsAnalysis_vectorMagnitude, 

    normalize, 
    vectorSubtract, 
    dot, 
    magnitude, 
    setMagnitude,

    test_xyzlatlon, 
    xyzFromLatLon_rad, 
    xyzToLatLon_rad, 
    xyzFromLatLon_deg, 
    xyzToLatLon_deg, 
    RAD2DEG, 
    DEG2RAD,
    deg_latlonVectorAt_to_xyzVector
};
