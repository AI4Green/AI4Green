import { format, parseISO } from "date-fns";
import { toZonedTime } from "date-fns-tz";
import { Base64 } from "js-base64";

import { getDateLocale } from "config/i18n";

/**
 * Try to get the value from an object,
 * but if the key is not present, return the key name instead of undefined
 *
 * @param {*} o Source Dictionary to look in
 * @param {*} k Key
 * @returns
 */
export const getValueOrKey = (o, k) => o[k] ?? k;

/**
 * Try to find the key of an object inside an object
 * by the value of a property on the inner object.
 *
 * @param {*} o
 * @param {*} prop
 * @param {*} value
 * @returns
 */
export const findKeyByPropertyValue = (o, prop, value) =>
  // _ should be universal for discarding!

  Object.entries(o).find(([k, v]) =>
    v[prop] != null ? v[prop] === value : k === value,
  )?.[0];

export const Base64UrlToUtf8 = (input) =>
  (!!input && Base64.decode(input)) || null;

export const Base64UrlToJson = (input) => JSON.parse(Base64UrlToUtf8(input));

export const Utf8ToBase64Url = (input) =>
  (!!input && Base64.encodeURI(input)) || "";

/**
 * Reduce an array of objects to a keyed dictionary of those objects
 * in the form
 * ```
 *     {
 *         [object[keyProp]]: object
 *         ...
 *     }
 * ```
 * @param {object[]} data a list of objects
 * @param {any} [keyProp] the property of each object to use as a unique dictionary key.
 *
 * defaults to `"id"`
 */
export const toDictionary = (data, keyProp = "id") =>
  // TODO: make this even more like C#'s `Enumerable.ToDictionary()`
  data.reduce((acc, datum) => {
    acc[datum[keyProp]] = datum;
    return acc;
  }, {});

export const getFormattedDate = (isoDateString, culture) => {
  const date = parseISO(isoDateString);
  return format(toZonedTime(date, "UTC"), "Pp (z)", {
    locale: getDateLocale(culture),
    timeZone: "UTC",
  });
};
