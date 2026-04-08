import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tsconfigPaths from "vite-tsconfig-paths";
import path from "path";

export default defineConfig({
  plugins: [react(), tsconfigPaths()],
  build: {
    outDir: path.resolve(__dirname, "../Webapp/sources/static/dist"),
    emptyOutDir: true,
    sourcemap: false,
    rollupOptions: {
      input: path.resolve(__dirname, "index.html"),
    },
  },
  server: {
    port: 8000,
    https: false, // uses mkcert-generated certificate
    proxy: {
      // Optional Python backend proxy. Uncomment if you want it
      // "/api": {
      //   target: "http://localhost:8000",
      //   changeOrigin: true,
      //   secure: false
      // }
    },
  },
});
