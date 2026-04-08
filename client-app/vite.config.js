import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tsconfigPaths from "vite-tsconfig-paths";
import path from "path";

export default defineConfig({
  plugins: [react(), tsconfigPaths()],
  base: "/spa/",
  server: {
    port: 8000,
    // proxy: {
    //   "/api": {
    //     target: "http://127.0.0.1:", // use env variable here for deployment
    //     changeOrigin: true,
    //     secure: false,
    //   },
    // },
  },
});
