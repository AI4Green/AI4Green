import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tsconfigPaths from "vite-tsconfig-paths";
import path from "path";

export default defineConfig({
  plugins: [react(), tsconfigPaths()],
  base: "/spa/",
  build: {
    outDir: path.resolve(__dirname, "../Webapp/sources/static/spa"),
    emptyOutDir: true,
    manifest: true,

    rollupOptions: {
      input: {
        app: path.resolve(__dirname, "index.html"),
        coshhForm: path.resolve(__dirname, "src/entries/coshh-form.jsx"),
      },
    },
  },
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
