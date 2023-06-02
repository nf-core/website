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
            iconUrl: '/src/assets/marker-icon-2x-green.png',
            shadowUrl: '/src/assets/marker-shadow.png',
            iconSize: [25, 41],
            iconAnchor: [12, 41],
            popupAnchor: [1, -34],
            shadowSize: [41, 41],
        });
        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        }).addTo(m);
        let latlng: [number, number][] = [];
        console.log('locations', locations);
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
                            `<img src="/src/assets/contributors/colour/${marker.image}" title="${marker.name}" class="contributor_map_logo"></img>`
                    );
                latlng.push(marker.location);
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
</script>

<div style="height:500px;width:90%" use:mapAction />

<style lang="scss">
    .contributor_map_logo {
        height: 45px;
        margin-top: 0.25rem;
    }
</style>
