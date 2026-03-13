<script lang="ts">
    import { tileLayer, marker, map, Icon } from "leaflet";
    import "leaflet-fullscreen";
    import { onMount } from "svelte";

    // Dynamically load CSS only when component mounts
    onMount(async () => {
        await import("leaflet/dist/leaflet.css");
        await import("leaflet-fullscreen/dist/leaflet.fullscreen.css");
    });

    interface Props {
        locations?: {
            location: [number, number];
            name: string;
            url: string;
            image?: string;
        }[];
    }

    let { locations = [] }: Props = $props();

    let m;
    function createMap(container) {
        m = map(container, {
            minZoom: 1.49,
            fullscreenControl: true,
            maxBounds: [
                [-90, -180],
                [90, 180],
            ],
            maxBoundsViscosity: 1.0,
        }).setView([20, 32.09], 1.49);

        let greenIcon = new Icon({
            iconUrl: "/images/marker-icon-2x-green.png",
            shadowUrl: "/images/marker-shadow.png",
            iconSize: [25, 41],
            iconAnchor: [12, 41],
            popupAnchor: [1, -34],
            shadowSize: [41, 41],
        });
        tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        }).addTo(m);

        locations.map(function (locationMarker) {
            if (locationMarker != null) {
                const image = locationMarker.image
                    ? `<img src="/images/contributors/colour/${locationMarker.image}" title="${locationMarker.name}" class="contributor_map_logo"></img>`
                    : "";
                marker(locationMarker.location, {
                    icon: greenIcon,
                })
                    .addTo(m)
                    .bindPopup(
                        `<h6><a href="${
                            locationMarker.url.startsWith("/events/")
                                ? locationMarker.url
                                : locationMarker.name.replaceAll("/[^a-z]+/", "-")
                        }">${locationMarker.name}</a></h6>${image}`,
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

<svelte:window onresize={resizeMap} />

<div class="map m-auto" use:mapAction></div>

<style lang="scss">
    .map {
        width: 100%;
        height: 460px;
        max-width: 800px;
    }
    @media (max-width: 767.98px) {
        // md-breakpoint
        .map {
            height: 300px;
        }
    }
    :global(.contributor_map_logo) {
        max-height: 5.2rem;
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
