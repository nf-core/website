<script lang="ts">
    import * as L from 'leaflet';
    import 'leaflet/dist/leaflet.css';

    export let locations: {
        location: [number, number];
        name: string;
        url: string;
        image: string;
    }[] = [];

    let map;
    function createMap(container) {
        let m = L.map(container, { minZoom: 2 }).setView([15.505, 10.09], 2);
        let greenIcon = new L.Icon({
            iconUrl: '/_astro/marker-icon-2x-green.png',
            shadowUrl: '/_astro/marker-shadow.png',
            iconSize: [25, 41],
            iconAnchor: [12, 41],
            popupAnchor: [1, -34],
            shadowSize: [41, 41],
        });
        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        }).addTo(m);
        locations.map(function (marker) {
            if (marker != null) {
                L.marker(marker.location, {
                    icon: greenIcon,
                })
                    .addTo(m)
                    .bindPopup(
                        '<h6><a href="#' +
                            marker.name.replaceAll('/[^a-z]+/', '-') +
                            '">' +
                            marker.name +
                            '</a></h6>' +
                            `<img src="/_astro/contributors/colour/${marker.image}" title="${marker.name}" class="contributor_map_logo"></img>`
                    );
            }
        });

        return m;
    }

    function mapAction(container) {
        map = createMap(container);
        return {
            destroy: () => {
                map.remove();
            },
        };
    }
    function resizeMap() {
        if (map) {
            map.invalidateSize();
        }
    }
</script>

<svelte:window on:resize={resizeMap} />

<div class="map m-auto" use:mapAction />

<style lang="scss">
    @import '@styles/_variables.scss';
    .map {
        height: 480px;
        width: 90%;
    }
    @include media-breakpoint-down(md) {
        .map {
            height: 350px;
            width: 100%;
        }
    }
    :global(.contributor_map_logo) {
        max-height: 5rem;
        margin-top: 0.25rem;
    }
    :global(.leaflet-popup-content-wrapper) {
        max-height: 15rem;
        :global(.leaflet-popup-content) {
            display: flex;
            flex-direction: column;
            justify-content: center;
        }
    }
</style>
