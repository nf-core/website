<script lang="ts">
    import { tileLayer, marker, map, Icon } from 'leaflet';
    import 'leaflet/dist/leaflet.css';

    export let locations: {
        location: [number, number];
        name: string;
        url: string;
        image?: string;
    }[] = [];

    let m;
    function createMap(container) {
        m = map(container, { minZoom: 2 }).setView([15.505, 10.09], 2);
        let greenIcon = new Icon({
            iconUrl: '/images/marker-icon-2x-green.png',
            shadowUrl: '/images/marker-shadow.png',
            iconSize: [25, 41],
            iconAnchor: [12, 41],
            popupAnchor: [1, -34],
            shadowSize: [41, 41],
        });
        tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        }).addTo(m);
        locations.map(function (locationMarker) {
            const image = locationMarker.image
                ? `<img src="/images/contributors/colour/${locationMarker.image}" title="${locationMarker.name}" class="contributor_map_logo"></img>`
                : '';
            if (locationMarker != null) {
                marker(locationMarker.location, {
                    icon: greenIcon,
                })
                    .addTo(m)
                    .bindPopup(
                        '<h6><a href="#' +
                            locationMarker.name.replaceAll('/[^a-z]+/', '-') +
                            '">' +
                            locationMarker.name +
                            '</a></h6>' +
                            image,
                    );
            }
        });

        return m;
    }

    function mapAction(container) {
        m = createMap(container);
        return {
            destroy: () => {
                m.remove();
            },
        };
    }
    function resizeMap() {
        if (m) {
            m.invalidateSize();
        }
    }
</script>

<svelte:window on:resize={resizeMap} />

<div class="map m-auto" use:mapAction />

<style lang="scss">
    @import '../styles/_variables.scss';
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
