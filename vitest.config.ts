import codspeedPlugin from "@codspeed/vitest-plugin";
import { defineConfig } from "vitest/config";

export default defineConfig({
    plugins: [codspeedPlugin()],
    test: {
        // Keep V8 background work (GC, compilation) on the main thread so the
        // larger markdown benchmarks are not dominated by thread synchronisation.
        execArgv: ["--single-threaded"],
        benchmark: {
            include: ["bench/**/*.bench.ts"],
        },
    },
});
