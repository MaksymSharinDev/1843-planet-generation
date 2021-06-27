const util = require('./util');

// TODO: rename to createSecondaryProjection
exports.createProjection_lines = (projection, lines_xyz, {min, range}, centerLat=0, centerLon=0, planetRadius=1) => {
    let result = [];

    for (let r = 0; r < lines_xyz.length; r++) {
        let [x, y, z] = lines_xyz[r];
        let [nx, ny, nz] = projection([x, y, z], centerLat, centerLon, planetRadius);
        x = 2 * ((nx-min)/range) - 1;
        y = 2 * ((ny-min)/range) - 1;
        z = 2 * ((nz-min)/range) - 1;
        result.push([nx, ny, nz]);
    }

    return {result};
}

// TODO: allow r_elevation to be passed in to make a 3D elevation map
exports.createProjection = (projection, r_xyz, centerLat=0, centerLon=0, planetRadius=1) => {
    let result = [];

    let min = 99999999;
    let max = -99999999;
    
    for (let r = 0; r < r_xyz.length/3; r++) {
        let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        let [nx, ny, nz] = projection([x, y, z], centerLat, centerLon, planetRadius);
        result.push(nx);
        result.push(ny);
        result.push(nz);
        
        min = Math.min(nx, min);
        max = Math.max(nx, max);
        min = Math.min(ny, min);
        max = Math.max(ny, max);
        min = Math.min(nz, min);
        max = Math.max(nz, max);
    }

    let range = max - min;

    for (let r = 0; r < r_xyz.length/3; r++) {
        let [nx, ny, nz] = result.slice(3 * r, 3 * r + 3);
        result[3*r+0] = 2 * ((nx-min)/range) - 1;
        result[3*r+1] = 2 * ((ny-min)/range) - 1;
        result[3*r+2] = 2 * ((nz-min)/range) - 1;
    }

    return {result, max, min, range};
}

// TODO: something fishy is happening here...
// longitude < 0 or > 180 is acting funny

exports.mercator = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_rad(x, y, z);

    let nx = planetRadius * (lon - centerLon);
    let ny = planetRadius * Math.log(Math.tan((Math.PI/4) + (lat/2)));
    return [nx, ny, 0];
}

exports.equirectangular = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_rad(x, y, z);

    let nx = planetRadius * (lon - centerLon);
    let ny = planetRadius * (lat - centerLat);
    
    return [nx, ny, 0];
}

exports.sinusoidal = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_rad(x, y, z);

    return [
        (lon - centerLon) * Math.cos(lat),
        lat, 
        0
    ];
}

exports.gall_stereographic = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_deg(x, y, z);

    lat = (Math.PI/180) * lat;
    lon = (Math.PI/180) * lon;

    return [
        planetRadius * lon * Math.SQRT1_2,
        planetRadius * (1 + Math.SQRT2 / 2) * Math.tan(lat / 2), 
        0
    ];
}

// credit for this projection and code: https://web.archive.org/web/20120407111351/http://www.shadedrelief.com/NE_proj/index.html
exports.natural_earth = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_deg(x, y, z);

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

    let nx = lon * (A0 + phi2 * (A1 + phi2 * (A2 + phi4 * phi2 * (A3 + phi2 * A4))));
    let ny = lat * (B0 + phi2 * (B1 + phi4 * (B2 + B3 * phi2 + B4 * phi4)));

    return [
        nx,
        ny,
        0
    ];
}