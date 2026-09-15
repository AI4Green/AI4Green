import { enGB } from "date-fns/locale";
import i18n from "i18next";
import LanguageDetector from "i18next-browser-languagedetector";
import Backend from "i18next-http-backend";
import { initReactI18next } from "react-i18next";

// we should preload all real langs
// as sometimes we want a string from NOT the current language
// e.g. for the language select UI
const preload = []; // TODO: AI4Green4Students doesn't use language select UI at this time

i18n
  .use(Backend)
  .use(LanguageDetector)
  .use(initReactI18next)
  .init({
    debug: import.meta.env.DEV,

    preload,

    backend: {
      loadPath: "/locales/{{lng}}/{{ns}}.json",
    },

    supportedLngs: [...preload, "dev"],
    nonExplicitSupportedLngs: true,

    defaultNS: "core",
    ns: [], // only load namespaces on demand

    // currently we only accommodate generic languages, no regional specialisations at this time
    // this setting simplifies that when using the http backend
    load: "languageOnly",

    interpolation: {
      escapeValue: false,
    },
  });

// we use a locale supporting date library
// but we need to opt in to locales
// so we currently maintain a manual lookup here
// based on the languages we support
export const getDateLocale = (culture) => {
  switch (culture) {
    default:
      return enGB;
  }
};

export default i18n;
