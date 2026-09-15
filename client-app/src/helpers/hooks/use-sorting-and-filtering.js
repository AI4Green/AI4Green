import { useEffect, useState } from "react";

/**
 * Gets the appropriate sort function for a given object property.
 *
 * @param {string} key The object property to sort by.
 * @param {boolean} asc Sort ascending or descending.
 * @param {object} sorters Custom sorting behaviours.
 *
 * It's possible to specify for a given sorting key an alternative object property to sort against,
 * or alternative sorting logic as a sorting predicate (asc) => (a, b) => {}
 *
 * example:
 * {
 *   active: { property: "activeInstanceId" },
 *   name: {sorter: asc => (a, b) => asc ? a.localeCompare(b) : b.localeCompare(a) }
 * }
 *
 */
const getPropertySorter = (key, asc, sorters = {}) => {
  const customSorter = sorters[key] ?? {};

  const { property, sorter } = customSorter;
  return ({ [property ?? key]: a }, { [property ?? key]: b }) => {
    if (sorter) {
      return sorter(asc)(a, b);
    } else {
      return asc ? a - b : b - a;
    }
  };
};

/**
 * Build a lookup list of object IDs and names
 * sorted by the specified object property,
 * ascending or descending.
 *
 * @param {object} input The source dictionary of objects.
 * @param {string} key The object property to sort by.
 * @param {boolean} asc Sort ascending or descending.
 */
const getSortedLookup = (input, key, asc, sorters) =>
  Object.keys(input)
    .map((id) => input[id])
    .sort(getPropertySorter(key, asc, sorters));

/**
 * Filters a list of objects with a `name` (or custom specified property name) property
 * where said property contains a match on the `filter` string.
 *
 * If no filter is specified, the full input list is returned.
 *
 * @param {object[]} input The source list of objects.
 * @param {string} filter The filter string to match against.
 * @param {string|function} filterer
 * If a string, is a property name used to match against filter
 *
 * If a function, is executed to perform custom filtering, with the signature (input, filter) => {}
 */
const getFilteredLookup = (input, filter, filterer = "name") =>
  !filter
    ? input
    : input.filter(
        typeof filterer === "string"
          ? ({ [filterer]: keyProperty }) =>
              new RegExp(filter, "i").test(keyProperty)
          : (value) => filterer(value, filter),
      );

const storageKeyPrefix = "ai4green4students.sorting";

export const useSortingAndFiltering = (
  sourceList,
  filterer,
  { initialSort, sorters, storageKey },
) => {
  const storageKeyFull = `${storageKeyPrefix}.${storageKey}`;

  const [sorting, setSorting] = useState(
    JSON.parse(localStorage.getItem(storageKeyFull)) || {
      key: initialSort?.key ?? "name",
      [initialSort.key ?? "name"]: initialSort?.asc,
    },
  );

  const [sortedList, setSortedList] = useState([]);
  const [filter, setFilter] = useState("");
  const [outputList, setOutputList] = useState([]);

  // update the sorted list appropriately
  useEffect(() => {
    setSortedList(
      getSortedLookup(sourceList, sorting.key, sorting[sorting.key], sorters),
    );
    // also save the new sort for next time
    localStorage.setItem(storageKeyFull, JSON.stringify(sorting));
  }, [sourceList, sorting, sorters, storageKeyFull]);

  // update the filtered list appropriately
  useEffect(() => {
    setOutputList(getFilteredLookup(sortedList, filter, filterer));
  }, [filter, filterer, sortedList]);

  // default sort handler
  const handleSort = (key) => {
    setSorting({
      ...sorting,
      key,
      [key]: sorting.key === key ? !sorting[key] : sorting[key],
    });
  };

  return {
    sorting,
    setSorting,
    onSort: handleSort,
    filter,
    setFilter,
    outputList,
  };
};
