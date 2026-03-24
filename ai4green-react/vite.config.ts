import { defineConfig } from "vite";
import react, { reactCompilerPreset } from "@vitejs/plugin-react";
import babel from "@rolldown/plugin-babel";

export default defineConfig({
  plugins: [react(), babel({ presets: [reactCompilerPreset()] })],
  resolve: {
    alias: {
      src: "/src",
    },
  },
  server: {
    proxy: {
      "/coshh_setup": "http://127.0.0.1:80",
    },
  },
});
