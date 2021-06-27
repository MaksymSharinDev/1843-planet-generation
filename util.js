
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
    return [Math.cos(latRad) * Math.sin(lonRad),
            Math.sin(latRad),
            Math.cos(latRad) * Math.cos(lonRad)];
}

function xyzToLatLon_rad(x, y, z) {
    // formulas from https://gis.stackexchange.com/a/120685
    let r = magnitude([x,y,z]);

    return [
        Math.asin(y / r),
        Math.atan2(x, z)
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

function test() {
    let latlon = [
        [0, 0],
        [10, 0],
        [20, 0],
        [-10, 0],
        [-20, 0],

        [90, 0],
        [-90, 0],
        [180, 0], // any lattitude greater than 180 or less than -180 is invalid, but these functions handle these cases very gracefully
        [-180, 0], // and interpret them generously to mean [0, 180]. beautiful :)
        
        [0, 10],
        [0, -10],
        [0, 360],
        [0, 370],
    ];

    let testU = ([lat, lon]) => {
        let [x, y, z] = xyzFromLatLon_deg(lat, lon);

        console.log([lat, lon] + " -> " + [x, y, z] + " -> " + xyzToLatLon_deg(x, y, z));
    }

    latlon.forEach(e => testU(e));
}

module.exports = { test, vectorSubtract, dot, magnitude, xyzFromLatLon_rad, xyzToLatLon_rad,  xyzFromLatLon_deg, xyzToLatLon_deg, RAD2DEG, DEG2RAD};
