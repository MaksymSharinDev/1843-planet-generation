/*
 * From http://www.redblobgames.com/x/1742-webgl-mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

/* Generate the biome colormap indexed by elevation -1:+1 and rainfall 0:1 */
exports.width = 64,
exports.height = 64;

function colormap_standard() {
    const pixels = new Uint8Array(exports.width * exports.height * 4);

    for (var y = 0, p = 0; y < exports.height; y++) {
        for (let x = 0; x < exports.width; x++) {
            let e = 2 * x / exports.width - 1,
                m = y / exports.height;
            
            let r, g, b;

            if (x === exports.width/2 - 1) {
                r = 48;
                g = 120;
                b = 160;
            } else
            if (x === exports.width/2 - 2) {
                r = 48;
                g = 100;
                b = 150;
            } else if (x === exports.width/2 - 3) {
                r = 48;
                g = 80;
                b = 140;
            } else
            if (e < 0.0) {
                r = 48 + 48*e;
                g = 64 + 64*e;
                b = 127 + 127*e;
            } else { // adapted from terrain-from-noise article
                m = m * (1-e); // higher elevation holds less moisture; TODO: should be based on slope, not elevation
                
                r = 210 - 100*m;
                g = 185 - 45*m;
                b = 139 - 45*m;
                r = 255 * e + r * (1-e),
                g = 255 * e + g * (1-e),
                b = 255 * e + b * (1-e);
            }

            pixels[p++] = r;
            pixels[p++] = g;
            pixels[p++] = b;
            pixels[p++] = 255;
        }
    }
    return pixels;
}

exports.colormap_standard = colormap_standard();

function colormap_temperature() {
    const pixels = new Uint8Array(exports.width * exports.height * 4);

    for (var y = 0, p = 0; y < exports.height; y++) {
        for (let x = 0; x < exports.width; x++) {
            let r, g, b;

            r = 255*y/exports.height;
            g = 0;
            b = 0;

            pixels[p++] = r;
            pixels[p++] = g;
            pixels[p++] = b;
            pixels[p++] = 255;
        }
    }
    return pixels;
}

exports.colormap_temperature = colormap_temperature();

function colormap_humidity() {
    const pixels = new Uint8Array(exports.width * exports.height * 4);

    for (var y = 0, p = 0; y < exports.height; y++) {
        for (let x = 0; x < exports.width; x++) {
            let r, g, b;

            r = 0;
            g = 0;
            b = 255*y/exports.height;

            pixels[p++] = r;
            pixels[p++] = g;
            pixels[p++] = b;
            pixels[p++] = 255;
        }
    }
    return pixels;
}

exports.colormap_humidity = colormap_humidity();

function colormap_cloudcover() {
    const pixels = new Uint8Array(exports.width * exports.height * 4);

    for (var y = 0, p = 0; y < exports.height; y++) {
        for (let x = 0; x < exports.width; x++) {
            let r, g, b, a;

            let prop = (y/exports.height);
            // desmos: \min\left(1,\frac{0.9}{1+e^{14\left(x-0.75\right)}}+0.12\right)
            let sigmoid = 0.8 / ( 1 + Math.exp(14 * (prop - 0.75)) ) + 0.22;
            r = g = b = 255 * Math.min(1, sigmoid);

            pixels[p++] = r;
            pixels[p++] = g;
            pixels[p++] = b;
            pixels[p++] = Math.min(255, 2.5*255*y/exports.height);
        }
    }
    return pixels;
}

exports.colormap_cloudcover = colormap_cloudcover();

function colormap_temperature_and_humidity() {
    const pixels = new Uint8Array(exports.width * exports.height * 4);

    for (var y = 0, p = 0; y < exports.height; y++) {
        for (let x = 0; x < exports.width; x++) {
            let r, g, b, a;

            let t_prop = x/exports.width;
            let h_prop = y/exports.width;
            
            r = 255 * (t_prop);
            g = 0;
            b = 255 * (h_prop);

            pixels[p++] = r;
            pixels[p++] = g;
            pixels[p++] = b;
            pixels[p++] = 255;
        }
    }
    return pixels;
}

exports.colormap_temperature_and_humidity = colormap_temperature_and_humidity();

function colormap_debug_solid_blue() {
    const pixels = new Uint8Array(exports.width * exports.height * 4);

    for (var y = 0, p = 0; y < exports.height; y++) {
        for (let x = 0; x < exports.width; x++) {
            let r, g, b, a;
            
            r = 0;
            g = 0;
            b = 255;

            pixels[p++] = r;
            pixels[p++] = g;
            pixels[p++] = b;
            pixels[p++] = 255;
        }
    }
    return pixels;
}

exports.colormap_debug_solid_blue = colormap_debug_solid_blue();