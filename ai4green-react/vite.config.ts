import { defineConfig } from "vite";
import react, { reactCompilerPreset } from "@vitejs/plugin-react";
import babel from "@rolldown/plugin-babel";
import path from "path";

export default defineConfig({
  plugins: [react(), babel({ presets: [reactCompilerPreset()] })],
  resolve: {
    alias: {
      src: path.resolve(__dirname, "src/"),
      components: path.resolve(__dirname, "src/components"),
      assets: path.resolve(__dirname, "src/assets"),
      api: path.resolve(__dirname, "src/api"),
      contexts: path.resolve(__dirname, "src/contexts"),
      pages: path.resolve(__dirname, "src/pages"),
      constants: path.resolve(__dirname, "src/constants"),
      layouts: path.resolve(__dirname, "src/layouts"),
      helpers: path.resolve(__dirname, "src/helpers"),
      config: path.resolve(__dirname, "src/config"),
      routes: path.resolve(__dirname, "src/routes"),
    },
  },
  server: {
    proxy: {
      "/coshh_setup": "http://127.0.0.1:80",
    },
  },
});
