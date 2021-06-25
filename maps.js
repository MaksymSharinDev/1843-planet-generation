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

function normalize(xyz) {
    let [x, y, z] = xyz;
    let mag = magnitude(xyz);

    return [x/mag, y/mag, z/mag];
}

function xyzToLatLon(xyz, planetRadius = 1) {
    let [x, y, z] = normalize(xyz);
    let south = false;
    if (z < 0) { z = -z; south = true; }

    let lat_deg = (180/Math.PI) * Math.acos(Math.abs(z) / planetRadius), 
        lon_deg = (180/Math.PI) * Math.atan2(y, x);

    if (south) lat_deg = 180-lat_deg;
    return [lat_deg, lon_deg];
}

exports.createProjection_lines = (projection, lines_xyz, {min, range}, centerLat=0, centerLon=0, planetRadius=1) => {
    let result = [];

    for (let r = 0; r < lines_xyz.length; r++) {
        let xyz = lines_xyz[r];
        let [x, y, z] = projection(xyz, centerLat, centerLon, planetRadius);
        x = 2 * ((x-min)/range) - 1;
        y = 2 * ((y-min)/range) - 1;
        z = 2 * ((z-min)/range) - 1;
        result.push([x, y, z]);
    }

    return {result};
}

// TODO: allow r_elevation to be passed in to make a 3D elevation map
exports.createProjection = (projection, r_xyz, centerLat=0, centerLon=0, planetRadius=1) => {
    let result = [];

    let min = 99999999;
    let max = -99999999;
    
    for (let r = 0; r < r_xyz.length/3; r++) {
        let xyz = r_xyz.slice(3 * r, 3 * r + 3);
        let [x, y, z] = projection(xyz, centerLat, centerLon, planetRadius);
        result.push(x);
        result.push(y);
        result.push(z);
        
        min = Math.min(x, min);
        max = Math.max(x, max);
        min = Math.min(y, min);
        max = Math.max(y, max);
        min = Math.min(z, min);
        max = Math.max(z, max);
    }

    let range = max - min;

    for (let r = 0; r < r_xyz.length/3; r++) {
        let [x, y, z] = result.slice(3 * r, 3 * r + 3);
        result[3*r+0] = 2 * ((x-min)/range) - 1;
        result[3*r+1] = 2 * ((y-min)/range) - 1;
        result[3*r+2] = 2 * ((z-min)/range) - 1;
    }

    return {result, max, min, range};
}

exports.mercator = (xyz, centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = xyzToLatLon(xyz, planetRadius);

    lat = (Math.PI/180) * lat;
    lon = (Math.PI/180) * lon;

    let x = planetRadius * (lon - centerLon);
    let y = planetRadius * Math.log(Math.tan((Math.PI/4) + (lat/2)));
    return [x, y, 0];
}

exports.equirectangular = (xyz, centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = xyzToLatLon(xyz, planetRadius);

    lat = (Math.PI/180) * lat;
    lon = (Math.PI/180) * lon;

    let x = planetRadius * (lon - centerLon);
    let y = planetRadius * (lat - centerLat);
    return [x, y, 0];
}

exports.sinusoidal = (xyz, centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = xyzToLatLon(xyz, planetRadius);

    lat = (Math.PI/180) * (lat+90);
    lon = (Math.PI/180) * lon;

    return [
        (lon - centerLon) * Math.cos(lat),
        lat, 
        0
    ];
}

exports.gall_stereographic = (xyz, centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = xyzToLatLon(xyz, planetRadius);

    lat = (Math.PI/180) * lat;
    lon = (Math.PI/180) * lon;

    return [
        planetRadius * lon * Math.SQRT1_2,
        planetRadius * (1 + Math.SQRT2 / 2) * Math.tan(lat / 2), 
        0
    ];
}

// credit for this projection and code: https://web.archive.org/web/20120407111351/http://www.shadedrelief.com/NE_proj/index.html
exports.natural_earth = (xyz, centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = xyzToLatLon(xyz, planetRadius);

    lat = (Math.PI/180) * lat;
    lon = (Math.PI/180) * lon;

    let l = (p) => 0.8707 - 0.1319

    const A0 = 0.8707;
    const A1 = -0.131979;
    const A2 = -0.013791;
    const A3 = 0.003971;
    const A4 = -0.001529;
    const B0 = 1.007226;
    const B1 = 0.015085;
    const B2 = -0.044475;
    const B3 = 0.028874;
    const B4 = -0.005916;
    const C0 = B0;
    const C1 = 3 * B1;
    const C2 = 7 * B2;
    const C3 = 9 * B3;
    const C4 = 11 * B4;
    const EPS = 1e-11;
    const MAX_Y = 0.8707 * 0.52 * Math.PI;


    const phi2 = lat * lat;
    const phi4 = phi2 * phi2;

    let x = lon * (A0 + phi2 * (A1 + phi2 * (A2 + phi4 * phi2 * (A3 + phi2 * A4))));
    let y = lat * (B0 + phi2 * (B1 + phi4 * (B2 + B3 * phi2 + B4 * phi4)));

    return [
        x,
        y,
        0
    ];
}