
const colormap = require('./colormap');
const maps = require('./maps');
const {vec3, mat4} = require('gl-matrix');

const textureKey = document.getElementById("textureKey").getContext('2d');
const draw_textureKey = ({texture}) => {
    var imgData = textureKey.createImageData(colormap.width, colormap.height);
    for (let i = 0; i < imgData.data.length; i++) { imgData.data[i] = texture[i]; }
    textureKey.putImageData(imgData, 0, 0, 0, 0, colormap.width, colormap.height);
};

const shaders = {
    renderPoints: regl => { 
        return {
            frag: `
                precision mediump float;
                void main() {
                    gl_FragColor = vec4(0, 0, 0, 1);
                }
            `,
                
            vert: `
                precision mediump float;
                uniform mat4 u_projection;
                uniform float u_pointsize;
                attribute vec3 a_xyz;
                void main() {
                    gl_Position = u_projection * vec4(a_xyz, 1);
                    gl_PointSize = gl_Position.z > 0.0? 0.0 : u_pointsize;
                }
            `,
        
            depth: {
                enable: false,
            },
            
            uniforms: {
                u_projection: regl.prop('u_projection'),
                u_pointsize: regl.prop('u_pointsize'),
            },
        
            primitive: 'points',
            count: regl.prop('count'),
            attributes: {
                a_xyz: regl.prop('a_xyz'),
            },
        }
    },
    renderLines: regl => { 
        return {
            frag: `
                precision mediump float;
                uniform vec4 u_multiply_rgba, u_add_rgba;
                varying vec4 v_rgba;
                void main() {
                    gl_FragColor = v_rgba * u_multiply_rgba + u_add_rgba;
                }
            `,
        
            vert: `
                precision mediump float;
                uniform mat4 u_projection;
                attribute vec3 a_xyz;
                attribute vec4 a_rgba;
                varying vec4 v_rgba;
                void main() {
                    vec4 pos = u_projection * vec4(a_xyz, 1);
                    v_rgba = (-2.0 * pos.z) * a_rgba;
                    gl_Position = pos;
                }
            `,
        
            depth: {
                enable: false,
            },
            
            uniforms: {
                u_projection: regl.prop('u_projection'),
                u_multiply_rgba: regl.prop('u_multiply_rgba'),
                u_add_rgba: regl.prop('u_add_rgba'),
            },
        
            // Interesting, this property keeps lines behind the planet from drawing
            blend: {
                enable: true,
                func: {src: 'one', dst: 'one minus src alpha'},
                equation: {
                    rgb: 'add',
                    alpha: 'add'
                },
                color: [0, 0, 0, 0],
            },
            primitive: 'lines',
            count: regl.prop('count'),
            attributes: {
                a_xyz: regl.prop('a_xyz'),
                a_rgba: regl.prop('a_rgba'),
            },
        }
    },
    renderTriangles_withElevation: regl => { 
        return {
            blend: {
                enable: true,
                func: {
                    srcRGB: 'src alpha',
                    srcAlpha: 1,
                    dstRGB: 'one minus src alpha',
                    dstAlpha: 1
                },
                equation: {
                    rgb: 'add',
                    alpha: 'add'
                },
                color: [0, 0, 0, 0]
            },
        
            frag: `
                precision mediump float;
                uniform sampler2D u_colormap;
                varying vec2 v_tm;
                void main() {
                    float e = v_tm.x > 0.0? 0.5 * (v_tm.x * v_tm.x + 1.0) : 0.5 * (v_tm.x + 1.0);
                    gl_FragColor = texture2D(u_colormap, vec2(e, v_tm.y));
                }
            `,
        
            vert: `
                precision mediump float;
                uniform mat4 u_projection;
                uniform float u_radius;
                attribute vec3 a_xyz;
                attribute vec2 a_tm;
                varying vec2 v_tm;
                void main() {
                    v_tm = a_tm;
                    gl_Position = u_projection * vec4(u_radius * a_xyz, 1);
                }
            `,
        
            uniforms: {
                u_colormap: regl.prop('u_colormap'),
                u_projection: regl.prop('u_projection'),
                u_radius: regl.prop('u_radius'),
            },
        
            count: regl.prop('count'),
            attributes: {
                a_xyz: regl.prop('a_xyz'),
                a_tm: regl.prop('a_tm'),
            },
        }
    },
    renderTriangles_linear: regl => { 
        return {
            blend: {
                enable: true,
                func: {
                    srcRGB: 'src alpha',
                    srcAlpha: 1,
                    dstRGB: 'one minus src alpha',
                    dstAlpha: 1
                },
                equation: {
                    rgb: 'add',
                    alpha: 'add'
                },
                color: [0, 0, 0, 0]
            },
        
            frag: `
                precision mediump float;
                uniform sampler2D u_colormap;
                varying vec2 v_tm;
                void main() {
                    gl_FragColor = texture2D(u_colormap, vec2(v_tm.x, v_tm.y));
                }
            `,
        
            vert: `
                precision mediump float;
                uniform mat4 u_projection;
                uniform float u_radius;
                attribute vec3 a_xyz;
                attribute vec2 a_tm;
                varying vec2 v_tm;
                void main() {
                    v_tm = a_tm;
                    gl_Position = u_projection * vec4(u_radius * a_xyz, 1);
                }
            `,
        
            uniforms: {
                u_colormap: regl.prop('u_colormap'),
                u_projection: regl.prop('u_projection'),
                u_radius: regl.prop('u_radius'),
            },
        
            count: regl.prop('count'),
            attributes: {
                a_xyz: regl.prop('a_xyz'),
                a_tm: regl.prop('a_tm'),
            },
        }
    },
    renderIndexedTriangles_withElevation: regl => {
        return {
            frag: `
                #extension GL_OES_standard_derivatives : enable
                
                precision mediump float;
                
                uniform sampler2D u_colormap;
                uniform vec2 u_light_angle;
                uniform float u_inverse_texture_size, u_slope, u_flat, u_c, u_d, u_outline_strength;
                
                varying vec2 v_tm;
                void main() {
                    float e = v_tm.x > 0.0? 0.5 * (v_tm.x * v_tm.x + 1.0) : 0.5 * (v_tm.x + 1.0);
                    float dedx = dFdx(v_tm.x);
                    float dedy = dFdy(v_tm.x);
                    vec3 slope_vector = normalize(vec3(dedy, dedx, u_d * 2.0 * u_inverse_texture_size));
                    vec3 light_vector = normalize(vec3(u_light_angle, mix(u_slope, u_flat, slope_vector.z)));
                    float light = u_c + max(0.0, dot(light_vector, slope_vector));
                    float outline = 1.0 + u_outline_strength * max(dedx,dedy);
                    gl_FragColor = vec4(texture2D(u_colormap, vec2(e, v_tm.y)).rgb * light / outline, 1);
                }
            `,
        
            vert: `
                precision mediump float;
                uniform mat4 u_projection;
                attribute vec3 a_xyz;
                attribute vec2 a_tm;
                varying vec2 v_tm;
                void main() {
                    v_tm = a_tm;
                    gl_Position = u_projection * vec4(a_xyz, 1);
                }
            `,
        
            uniforms: {
                u_colormap: regl.prop('u_colormap'),
                u_projection: regl.prop('u_projection'),
                u_light_angle: [Math.cos(Math.PI/3), Math.sin(Math.PI/3)],
                u_inverse_texture_size: 1.0 / 2048,
                u_d: 60,
                u_c: 0.15,
                u_slope: 6,
                u_flat: 2.5,
                u_outline_strength: 5,
            },
        
            elements: regl.prop('elements'),
            attributes: {
                a_xyz: regl.prop('a_xyz'),
                a_tm: regl.prop('a_tm'),
            },
        }
    }
};

const r_color_functions = {
    surface: (map) => (r) => {
        let m = map.r_moisture[r];
        let e = map.r_elevation[r];
        return [e, m];
    },

    temperature: (map) => (r) => {
        let h = map.r_humidity[r];
        let t = map.r_temperature[r];
        return [0, t];
    },

    humidity: (map) => (r) => {
        let h = map.r_humidity[r];
        let t = map.r_temperature[r];
        return [0, h];
    },

    cloudcover: (map) => (r) => {
        let h = map.r_clouds[r];
        return [0, h];
    },

    temperature_and_humidity: (map) => (r) => {
        return [map.r_temperature[r], map.r_humidity[r]];
    },
};

// TODO: add a colormap that has muted surface colors so that lines are easier to see
// start with one that makes all colors much darker
const colormaps = {
    surface: {
        width: colormap.width,
        height: colormap.height,
        data: colormap.colormap_standard,
        wrapS: 'clamp',
        wrapT: 'clamp'
    },
    temperature: {
        width: colormap.width,
        height: colormap.height,
        data: colormap.colormap_temperature,
        wrapS: 'clamp',
        wrapT: 'clamp'
    },
    humidity: {
        width: colormap.width,
        height: colormap.height,
        data: colormap.colormap_humidity,
        wrapS: 'clamp',
        wrapT: 'clamp'
    },
    cloudcover: {
        width: colormap.width,
        height: colormap.height,
        data: colormap.colormap_cloudcover,
        wrapS: 'clamp',
        wrapT: 'clamp'
    },
    temperature_and_humidity: {
        width: colormap.width,
        height: colormap.height,
        data: colormap.colormap_temperature_and_humidity,
        wrapS: 'clamp',
        wrapT: 'clamp'
    }
};
const renderEnvironments = {
    // textureKey: {
    //     draw: draw_textureKey,
    //     drawByName: (name) => { draw_textureKey(colormaps[name].data); },
    // },

    globe: {
        regl: require('regl')({
            canvas: "#output",
            extensions: ['OES_element_index_uint', 'OES_standard_derivatives']
        }),
        textures: {},
        shaders: {},
        r_color_functions: {},
        
        state: {
            lastProjection: { projectionFunction: null, projectionData: null },
        },
    },
    map: {
        regl: require('regl')({
            canvas: "#map",
            extensions: ['OES_element_index_uint', 'OES_standard_derivatives']
        }),
        textures: {},
        shaders: {},
        r_color_functions: {},
        
        state: {
            lastProjection: { projectionFunction: null, projectionData: null },
        },
    },
}

function initRenderEnvironments() {
    for ([envName, environment] of Object.entries(renderEnvironments)) {
        if (envName === "textureKey") continue; // special case

        for ([textureName, texture] of Object.entries(colormaps)) {
            environment.textures[textureName] = environment.regl.texture(texture);
        }
        for ([shaderName, shader] of Object.entries(shaders)) {
            environment.shaders[shaderName] = environment.regl(shader(environment.regl));
        }
        for ([layerName, r_color_function] of Object.entries(r_color_functions)) {
            environment.r_color_functions[layerName] = r_color_function;
        }
    }
}
initRenderEnvironments();

const drawFunctions = {
    drawLattitudeLine: (u_projection, latDeg, renderFunction, xyzSecondaryProjection, projectionData, lonStepDeg = 2, color = [1, 0, 0, 0]) => {
        if (Math.abs(latDeg) >= 90) return;
    
        let line_xyz = [], line_rgba = [];
    
        let latRad     = latDeg     / 180.0 * Math.PI,
            lonStepRad = lonStepDeg / 180.0 * Math.PI;
        let lonRad = 0;
        
        let lastPoint = [Math.cos(latRad) * Math.cos(lonRad),
                         Math.cos(latRad) * Math.sin(lonRad),
                         Math.sin(latRad)];
    
        for (let lonRad = lonStepRad; lonRad <= 2*Math.PI; lonRad += lonStepRad) {
            let nextPoint =  [Math.cos(latRad) * Math.cos(lonRad),
                              Math.cos(latRad) * Math.sin(lonRad),
                              Math.sin(latRad)];
    
            line_xyz.push(lastPoint, nextPoint);
            line_rgba.push(color, color);
    
            lastPoint = nextPoint;
        }
    
        renderFunction({
            u_projection,
            u_multiply_rgba: [1, 1, 1, 1],
            u_add_rgba: [0, 0, 0, 0],
            a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result,
            a_rgba: line_rgba,
            count: line_xyz.length,
        });
    },

    drawLongitudeLine: (u_projection, lonDeg, renderFunction, xyzSecondaryProjection, projectionData, latStepDeg = 2, color = [1, 0, 0, 0]) => {
        if (Math.abs(lonDeg) >= 90) return;

        let line_xyz = [], line_rgba = [];

        let lonRad     = lonDeg     / 180.0 * Math.PI,
            latStepRad = latStepDeg / 180.0 * Math.PI;
        let latRad = 0.5*Math.PI;
        
        let lastPoint = [Math.cos(latRad) * Math.cos(lonRad),
                        Math.cos(latRad) * Math.sin(lonRad),
                        Math.sin(latRad)];

        for (let latRad = 0.5*Math.PI + latStepRad; latRad <= 1.5*Math.PI; latRad += latStepRad) {
            let nextPoint =  [Math.cos(latRad) * Math.cos(lonRad),
                            Math.cos(latRad) * Math.sin(lonRad),
                            Math.sin(latRad)];

            line_xyz.push(lastPoint, nextPoint);
            line_rgba.push(color, color);

            lastPoint = nextPoint;
        }

        renderFunction({
            u_projection,
            u_multiply_rgba: [1, 1, 1, 1],
            u_add_rgba: [0, 0, 0, 0],
            a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result,
            a_rgba: line_rgba,
            count: line_xyz.length,
        });
    },

    drawRivers: (u_projection, mesh, {t_xyz, s_flow}, drawFuncition, xyzSecondaryProjection, mapProjectionResults) => {
        let line_xyz = [], line_rgba = [];
    
        for (let s = 0; s < mesh.numSides; s++) {
            if (s_flow[s] > 1) {
                let flow = 0.1 * Math.sqrt(s_flow[s]);
                let inner_t = mesh.s_inner_t(s),
                    outer_t = mesh.s_outer_t(s);
                line_xyz.push(t_xyz.slice(3 * inner_t, 3 * inner_t + 3),
                              t_xyz.slice(3 * outer_t, 3 * outer_t + 3));
                if (flow > 1) flow = 1;
                let rgba_premultiplied = [0.2 * flow, 0.5 * flow, 0.7 * flow, flow];
                line_rgba.push(rgba_premultiplied, rgba_premultiplied);
            }
        }
    
        drawFuncition({
            u_projection,
            u_multiply_rgba: [1, 1, 1, 1],
            u_add_rgba: [0, 0, 0, 0],
            a_xyz: xyzSecondaryProjection(line_xyz, mapProjectionResults).result,
            a_rgba: line_rgba,
            count: line_xyz.length,
        });
    },

    drawRegionBoundaries: (mesh, map, u_projection, regionsHaveBoundary, renderLines, xyzSecondaryProjection, projectionData, color) => {
        let {t_xyz} = map;
        let line_xyz = [], line_rgba = [];
        for (let s = 0; s < mesh.numSides; s++) {
            let begin_r = mesh.s_begin_r(s),
                end_r = mesh.s_end_r(s);
            // if (r_plate[begin_r] !== r_plate[end_r]) {
            //if ((r_elevation[begin_r] <= 0) !== (r_elevation[end_r] <= 0)) {
            if (regionsHaveBoundary(map, begin_r, end_r)) {
                let inner_t = mesh.s_inner_t(s),
                    outer_t = mesh.s_outer_t(s);
                line_xyz.push(t_xyz.slice(3 * inner_t, 3 * inner_t + 3),
                                t_xyz.slice(3 * outer_t, 3 * outer_t + 3));
                line_rgba.push(color, color);
            }
        }
        renderLines({
            u_projection,
            u_multiply_rgba: [1, 1, 1, 1],
            u_add_rgba: [0, 0, 0, 0],
            a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result, 
            a_rgba: line_rgba,
            count: line_xyz.length,
        });
    }
};

let defaultOptions = {
    mesh: null,
    map: null,
    N: null,

    procession: 0,
    tilt: 0,
    rotation: 0,
    layer: "surface", // "temperature", "surfaceonly", ...
    generateVoronoiGeometry: null,

    cloud_height: 0,

    draw_axis: true,
    draw_plateVectors: false,
    draw_plateBoundaries: false,
    draw_landBoundaries: false,
    draw_elevationLines: false,
    draw_windVectors: false,
    draw_normalVectors: false,
    draw_equator: false,
    draw_extraLat: false,
    draw_primeMeridian: false,
    draw_points: false,

    forceRedraw: false,
    surfaceMode: "centroid",
    projection: null,
};

function extractOptions(options, envName) {
    let specificOptions = options[envName] || {};
    let genericOptions = options.all || {};

    let retval = {};
    for([optionName, defaultValue] of Object.entries(defaultOptions)) {
        let specific = specificOptions[optionName];
        let generic = genericOptions[optionName];
        let previous = previousOptions[envName] ? previousOptions[envName][optionName] : undefined;
        let defaultVal = defaultOptions[optionName];

        if (specific !== undefined) retval[optionName] = specific;
        else if (generic !== undefined) retval[optionName] = generic;
        else if (previous !== undefined) retval[optionName] = previous;
        else if (defaultVal !== undefined) retval[optionName] = defaultVal;
    }

    return retval;
}

let previousOptions = {}
function draw(options) {
    for ([envName, environment] of Object.entries(renderEnvironments)) {
        let envOptions = extractOptions(options, envName);

        if (!envOptions.forceRedraw && JSON.stringify(envOptions) === JSON.stringify(previousOptions[envName])) continue;
        if (!envOptions.map || !envOptions.mesh) continue;

        previousOptions[envName] = envOptions;
        //
        //
        // Setup
        //
        //

        let {mesh, map, N} = envOptions;

        let u_projection = mat4.create();
        mat4.scale(u_projection, u_projection, [1, 1, 0.5, 1]); // avoid clipping
        
        mat4.rotate(u_projection, u_projection, envOptions.procession, [1, 0, 0]);
        mat4.rotate(u_projection, u_projection, envOptions.tilt, [0, 1, 0]);
        mat4.rotate(u_projection, u_projection, envOptions.rotation, [0, 0, 1]);

        let texture = envOptions.layer === "surfaceonly" 
            ? environment.textures["surface"] 
            : environment.textures[envOptions.layer]; 
        
        let r_color_fn = envOptions.layer === "surfaceonly" 
            ? environment.r_color_functions.surface(map)
            : environment.r_color_functions[envOptions.layer](map);

        // TODO: add support in options for centerLat and centerLon
        let xyzProjection = envOptions.projection 
            ? (xyz) => maps.createProjection(envOptions.projection, xyz)
            : (xyz) => { return { result: xyz }; };
        
        let xyzSecondaryProjection = envOptions.projection
            ? (xyz, projectionData) => { return maps.createProjection_lines(envOptions.projection, xyz, projectionData); }
            : (xyz, projectionData) => { return { result: xyz }; };

        //
        //
        // Start drawing
        //
        //


        draw_textureKey({
            texture: envOptions.layer === 'surfaceonly' ? colormaps.surface : colormaps[envOptions.layer]
        });

        // TODO: support for quads mode
        let projectionData;
        if (true || drawMode === 'centroid') {
            let triangleGeometry = envOptions.generateVoronoiGeometry(mesh, map, r_color_fn);
  
            // Handle caching the projection
            projectionData = !environment.state.lastProjection.projectionData || environment.state.lastProjection.projectionFunction !== envOptions.projection
            ? xyzProjection(triangleGeometry.xyz)
            : environment.state.lastProjection.projectionData;
            
            environment.state.lastProjection.projectionData = projectionData;
            environment.state.lastProjection.projectionFunction = envOptions.projectionFunction;
            // end caching

            let shader = envOptions.layer === "surface" || envOptions.layer === "surfaceonly"
                ? environment.shaders.renderTriangles_withElevation
                : environment.shaders.renderTriangles_linear;

            shader({
                u_projection,
                u_colormap: texture,
                u_radius: 1,
                a_xyz: projectionData.result,
                a_tm: triangleGeometry.tm,
                count: triangleGeometry.xyz.length / 3,
            });

            drawFunctions.drawRivers(u_projection, mesh, map, environment.shaders.renderLines, xyzSecondaryProjection, projectionData);

            if (envOptions.layer === 'surface') {
                // recalculating triangleGeometry to get new triangleGeometry.tm specific to clouds
                
                // TODO: figure out why this code causes a crash
                // let cloudTriangleGeometry = envOptions.generateVoronoiGeometry(mesh, map, environment.r_color_functions.cloudcover);

                // environment.shaders.renderTriangles_linear({
                //     u_projection,
                //     u_colormap: environment.textures.cloudcover,
                //     u_radius: 1+envOptions.cloud_height,
                //     a_xyz: projectionData.result,
                //     a_tm: cloudTriangleGeometry.tm,
                //     count: cloudTriangleGeometry.xyz.length / 3,
                // });
            }
        } else if (drawMode === 'quads') {
            // renderIndexedTriangles({
            //     u_projection,
            //     a_xyz: quadGeometry.xyz,
            //     a_tm: quadGeometry.tm,
            //     elements: quadGeometry.I,
            // });

            // drawRivers(u_projection, mesh, map);

            // // recalculating triangleGeometry to get new triangleGeometry.tm specific to clouds
            // quadGeometry = generateQuadGeometry(mesh, map, environment.r_color_functions.cloudcover);
            // renderTriangles({
            //     u_projection,
            //     u_colormap: environment.textures.cloudcover,
            //     u_radius: 1.01,
            //     a_xyz: projectionData.result,
            //     a_tm: triangleGeometry.tm,
            //     count: triangleGeometry.xyz.length / 3,
            // });
            
            // draw_textureKey({
            //     texture: colormaps.surface
            // });
        }

        if (envOptions.draw_plateVectors) {
            let {r_xyz, r_plate, plate_vec} = map;
            let line_xyz = [], line_rgba = [];
        
            for (let r = 0; r < mesh.numRegions; r++) {
                line_xyz.push(r_xyz.slice(3 * r, 3 * r + 3));
                line_rgba.push([1, 0, 0, 1]);
                line_xyz.push(vec3.add([], r_xyz.slice(3 * r, 3 * r + 3),
                                    vec3.scale([], plate_vec[r_plate[r]], 2 / Math.sqrt(N))));
                line_rgba.push([1, 1, 1, 1]);
            }
        
            environment.shaders.renderLines({
                u_projection,
                u_multiply_rgba: [1, 1, 1, 1],
                u_add_rgba: [0, 0, 0, 0],
                a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result, 
                a_rgba: line_rgba,
                count: line_xyz.length,
            });
        }
        
        if (envOptions.draw_elevationLines) { 
            const elevationDiff = 0.3;
            drawFunctions.drawRegionBoundaries(
                mesh, 
                map, 
                u_projection, 
                ({r_elevation}, begin_r, end_r) => Math.floor(r_elevation[begin_r]/elevationDiff) !== Math.floor(r_elevation[end_r]/elevationDiff),
                environment.shaders.renderLines, 
                xyzSecondaryProjection, 
                projectionData, 
                [0.2, 0.7, 0.3, 1]
            );
        }

        if (envOptions.draw_landBoundaries || envOptions.draw_elevationLines) {
            drawFunctions.drawRegionBoundaries(
                mesh, 
                map, 
                u_projection, 
                ({r_elevation}, begin_r, end_r) => (r_elevation[begin_r] <= 0) !== (r_elevation[end_r] <= 0),
                environment.shaders.renderLines, 
                xyzSecondaryProjection, 
                projectionData, 
                [0, 1, 1, 1]
            );
        }

        // TODO: colorcode plate boundaries on sliding scale:
        // red: full convergent ----- white: neither ----- green: full divergent
        if (envOptions.draw_plateBoundaries) {
            drawFunctions.drawRegionBoundaries(
                mesh, 
                map, 
                u_projection, 
                ({r_plate}, begin_r, end_r) => r_plate[begin_r] !== r_plate[end_r], 
                environment.shaders.renderLines, 
                xyzSecondaryProjection, 
                projectionData, 
                [1, 0, 0, 1]
            );
        }

        if (envOptions.draw_windVectors) {
            let {r_xyz, r_wind} = map;
            let line_xyz = [], line_rgba = [];
        
            for (let r = 0; r < mesh.numRegions; r++) {
                line_xyz.push(r_xyz.slice(3 * r, 3 * r + 3));
                line_rgba.push([0.3, 0.3, 1, 1]);
                line_xyz.push(vec3.add([], r_xyz.slice(3 * r, 3 * r + 3),
                                    vec3.scale([], r_wind[r], 100 * 2 / Math.sqrt(N))));
                line_rgba.push([1, 1, 1, 1]);
            }
        
            environment.shaders.renderLines({
                u_projection,
                u_multiply_rgba: [1, 1, 1, 1],
                u_add_rgba: [0, 0, 0, 0],
                a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result,
                a_rgba: line_rgba,
                count: line_xyz.length,
            });
        }
        if (envOptions.draw_normalVectors) {
            let {r_xyz} = map;
            let line_xyz = [], line_rgba = [];
        
            for (let r = 0; r < mesh.numRegions; r++) {
                let [x, y, z] = r_xyz.slice(3 * r, 3 * r + 3);
        
                line_xyz.push([x, y, z]);
                line_rgba.push([0, 0, 0, 1]);
                line_xyz.push(vec3.add([], [x, y, z],
                                       vec3.scale([], [2*x, 2*y, 2*z], 2 / Math.sqrt(N))));
                line_rgba.push([0, 1, 0, 1]);
            }
        
            environment.shaders.renderLines({
                u_projection,
                u_multiply_rgba: [1, 1, 1, 1],
                u_add_rgba: [0, 0, 0, 0],
                a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result,
                a_rgba: line_rgba,
                count: line_xyz.length,
            });  
        }
    
        if (envOptions.draw_axis) {
            let line_xyz = [], line_rgba = [];

            let radius = 1;
        
            line_xyz.push([0,0,radius], [0,0,1.5*radius]);
            line_rgba.push([0,0,0,1], [1,0,0,1]);
        
            line_xyz.push([0,0,-radius], [0,0,-1.5*radius]);
            line_rgba.push([0,0,0,1], [1,0,0,1]);
                
            environment.shaders.renderLines({
                u_projection,
                u_multiply_rgba: [1, 1, 1, 1],
                u_add_rgba: [0, 0, 0, 0],
                a_xyz: xyzSecondaryProjection(line_xyz, projectionData).result,
                a_rgba: line_rgba,
                count: line_xyz.length,
            });
        }
        if (envOptions.draw_equator || envOptions.draw_extraLat) {
            drawFunctions.drawLattitudeLine(u_projection,  0, environment.shaders.renderLines, xyzSecondaryProjection, projectionData, 10, [0.6, 0, 0, 1]);
        }
        if (envOptions.draw_extraLat) {
            drawFunctions.drawLattitudeLine(u_projection,  30, environment.shaders.renderLines, xyzSecondaryProjection, projectionData, 10, [0.2, 0, 0, 1]);
            drawFunctions.drawLattitudeLine(u_projection, -30, environment.shaders.renderLines, xyzSecondaryProjection, projectionData, 10, [0.2, 0, 0, 1]);
            
            drawFunctions.drawLattitudeLine(u_projection,  60, environment.shaders.renderLines, xyzSecondaryProjection, projectionData, 10, [0.2, 0, 0, 1]);
            drawFunctions.drawLattitudeLine(u_projection, -60, environment.shaders.renderLines, xyzSecondaryProjection, projectionData, 10, [0.2, 0, 0, 1]);
        }
        if (envOptions.draw_primeMeridian) {
            drawFunctions.drawLongitudeLine(u_projection, 0, environment.shaders.renderLines, xyzSecondaryProjection, projectionData, 10, [0.6, 0, 0, 1]);
        }

        if (envOptions.draw_points) {
            environment.shaders.renderPoints({
                u_projection,
                u_pointsize: 0.1 + 100 / Math.sqrt(N),
                a_xyz: projectionData.result,
                count: mesh.numRegions,
            });
        }

        // TODO: test points, start with some xyz points that lie on the surface of the sphere
        // and convert them to lat, lon
        // draw both the lat lon lines and the xyz points separately
        // then do the reverse (start with lat, lon)
        // and then export the xyz to lat lon functions
    }
}


module.exports = {renderEnvironments, draw};