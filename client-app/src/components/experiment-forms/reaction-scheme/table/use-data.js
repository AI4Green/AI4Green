import { useField } from "formik";
import { useEffect, useMemo, useState } from "react";

export const useReactionTable = (name, ketcherData) => {
  const [field, , helpers] = useField(name);
  const hasExistingTableData = field.value?.length >= 1;

  const { initialTableData } = useInitialReactionTableData(ketcherData);

  const initialData = hasExistingTableData ? field.value : initialTableData;

  const [tableData, setTableData] = useState(initialData);

  useEffect(() => {
    helpers.setValue(tableData);
  }, [helpers, tableData]);

  return { tableData, setTableData };
};

const useInitialReactionTableData = (reactionData) => {
  const initialTableData = useMemo(
    () =>
      reactionData?.map((data) => ({
        substanceType: data?.substanceType,
        substancesUsed: data?.name,
        molWeight: data?.molecularWeight,
        density: data?.density,
        hazards: data?.hazards,
        mass: { value: 0, unit: "" },
      })),
    [reactionData],
  );

  return { initialTableData: initialTableData ?? [] };
};
