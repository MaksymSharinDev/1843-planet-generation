const util = require('./util');

exports.createProjection_lines = (projection, lines_xyz, {min, range}, centerLat=0, centerLon=0, planetRadius=1) => {
    let result = [];

    for (let r = 0; r < lines_xyz.length; r++) {
        let [x, y, z] = lines_xyz[r];
        let [nx, ny, nz] = projection([x, y, z], centerLat, centerLon, planetRadius);
        nx = 2 * ((nx-min)/range) - 1;
        ny = 2 * ((ny-min)/range) - 1;
        nz = 2 * ((nz-min)/range) - 1;
        result.push([nx, ny, nz]);
    }

    // remove lines that wrap around a seam in the map
    // not perfect, but that's ok
    for (let r = 0; r < result.length; r += 2) {
        let xyz = result[r];
        let xyz2 = result[r+1];

        if (util.magnitude(xyz, xyz2) > 1) {
            result[r+1] = [...xyz];
        }
    }

    return {result};
}

// TODO: allow r_elevation to be passed in to make a 3D elevation map
exports.createProjection = (projection, r_xyz, centerLat=0, centerLon=0, planetRadius=1) => {
    let result = [];

    let min = 99999999;
    let max = -99999999;
    
    // let test = [];

    let badTriangles = []; // foreach (x) => x and the next 8 elements are bad data (ie, 3 xyz triples starting at idx x)
    let lastxyz;
    let lastsrcxyz;
    let count = 0;
    // let output = [];
    
    for (let r = 0; r < r_xyz.length/3; r++) {
        let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        let [nx, ny, nz] = projection([x, y, z], centerLat, centerLon, planetRadius);
        
        
        min = Math.min(nx, min);
        max = Math.max(nx, max);
        min = Math.min(ny, min);
        max = Math.max(ny, max);
        min = Math.min(nz, min);
        max = Math.max(nz, max);

        // if (count > 0 && util.magnitude(lastxyz, [nx, ny, nz]) > 1) {
        // if (count > 0){
            // let lon1 = util.xyzToLatLon_deg(x, y, z)[1];
            // let lon2 = util.xyzToLatLon_deg(...lastxyz)[1];
            // if (count > 0 && Math.abs(lon1 - lon2) > 179) {
            //     // badTriangles.push(r-count);
            //     nx = lastxyz[0];
            //     ny = lastxyz[1];
            //     nz = lastxyz[2];
            //     console.log("bad triangle detected");
            //     console.log(lon1 + " ~ " + lon2)
            //     // output.push([x, y, z] + " ; " + lastsrcxyz + " -> " + util.xyzToLatLon_deg(x, y, z) + " ; " + util.xyzToLatLon_deg(...lastsrcxyz) + " -> " + [nx, ny, nz] + " ; " + lastxyz);
            // }
        // }
        // lastxyz = [nx, ny, nz];
        // lastsrcxyz = [x, y, z];
        // count = (count+1) % 3;



        result.push(nx);
        result.push(ny);
        result.push(nz);



        // if (-0.1 < nx && nx < 0.1) {
        //     xdata.push(x);
        //     ydata.push(y);
        //     zdata.push(z);
        // }

        // test[r] = [x, y, z] + " -> ";
    }

    // console.log(output);

    // console.log(badTriangles.map((i) => {
    //     let xyz1 = r_xyz.slice(3 * 3 * i + 0, 3 * 3 * i + 3);
    //     let xyz2 = r_xyz.slice(3 * 3 * i + 3, 3 * 3 * i + 6);
    //     let xyz3 = r_xyz.slice(3 * 3 * i + 6, 3 * 3 * i + 9);
        
    //     return [xyz1, xyz2, xyz3];
    // }));

    // util.statsAnalysis(xdata);
    // util.statsAnalysis(ydata);
    // util.statsAnalysis(zdata);

    // console.log(util.xyzToLatLon_deg(...[ 0.7104270818768333, -0.447447726674183, -0.542612513582856 ]));
    // console.log(util.xyzToLatLon_deg(...[ 0.7206956934886738, -0.42247528271210677, -0.5491164106541863 ]));
    // console.log(util.xyzToLatLon_deg(...[ 0.7033963847193084, -0.44069237235793246, -0.5576950411374877 ]));
    
    

    let range = max - min;

    for (let r = 0; r < r_xyz.length/3; r++) {
        let [nx, ny, nz] = result.slice(3 * r, 3 * r + 3);
        result[3*r+0] = 2 * ((nx-min)/range) - 1;
        result[3*r+1] = 2 * ((ny-min)/range) - 1;
        result[3*r+2] = -1;//2 * ((nz-min)/range) - 1;

        // test[r] += [result[3*r], result[3*r+1], result[3*r+2]];
    }

    // console.log(test);
    // console.log(util.xyzToLatLon_deg(...util.xyzFromLatLon_deg(0, -10)))

    return {result, max, min, range};
}

// TODO: something fishy is happening here...
// longitude < 0 or > 180 is acting funny

exports.mercator = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_rad(x, y, z);

    lon = lon - centerLon;
    // lon = (lon + (2*Math.PI)) % (2*Math.PI); 

    // lat = (lat + (0.5*Math.PI)) % (2*Math.PI);

    let nx = planetRadius * lon;
    let ny = planetRadius * Math.log(Math.tan((Math.PI/4) + (lat/2)));
    return [nx, ny, 0];
}

let count = 0;
exports.equirectangular = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_deg(x, y, z);

    // lon = lon > Math.PI / 2 ? lon + Math.PI : lon;
    // lon = lon < -Math.PI / 2 ? lon - Math.PI : lon;
    lon = lon - centerLon*util.RAD2DEG

    // lon = (lon+360) % 360;
    // lat = (lat+90) % 360;

    let nx = planetRadius * (lon);
    let ny = planetRadius * (lat);// - centerLat*util.RAD2DEG);
    
    // if (lon>358 /*&& count++ < 5*/) {
    //     // console.log([x, y, z] + " | " +lon + " | " + nx);
        
    //     ny = 500
    //     // nx = 180
    //     // if () crash
    // }

    return [nx, ny, 0];
}

exports.sinusoidal = ([x, y, z], centerLat = 0, centerLon = 0, planetRadius = 1) => {
    let [lat, lon] = util.xyzToLatLon_rad(x, y, z);

    // the below line proves that the issue is with 
    // xyzToLatLon_rad - it has issues producing negative longitudes
    // probably something to do with creating 0 instead of -180
    lon = -lon;

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