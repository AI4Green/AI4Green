import { createSystem, defaultConfig } from "@chakra-ui/react";

export const theme = createSystem(defaultConfig, {
  theme: {
    tokens: {
      fonts: {
        heading: { value: "Public Sans Variable, sans-serif" },
        body: { value: "Public Sans Variable, sans-serif" },
      },

      colors: {
        blue: {
          500: { value: "#2D78BE" },
        },
        green: {
          500: { value: "#2E8456" },
        },
        brand: {
          500: { value: "#2E7E45" },
        },
      },

      fontSizes: {
        sm: { value: "0.875rem" },
        md: { value: "1rem" },
        lg: { value: "1.125rem" },
      },

      fontWeights: {
        normal: { value: 400 },
        medium: { value: 500 },
        bold: { value: 700 },
      },

      shadows: {
        callout: {
          value: "0 2px 10px 0 rgba(0,0,0,.132)",
        },
      },
    },
  },
});
