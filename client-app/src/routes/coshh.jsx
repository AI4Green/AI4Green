import { Routes, Route } from "react-router-dom";
import {
  RoutedCoshhCreateModal,
  EmbeddedCoshh,
} from "components/coshh-forms/form.jsx";
import { CoshhInstanceList } from "pages/coshh-forms";

export const COSHH = () => {
  return (
    <Routes>
      <Route path="my-coshh-forms" element={<CoshhInstanceList />} />

      <Route path="new" element={<RoutedCoshhCreateModal />} />

      <Route path="form/:formId/edit" element={<EmbeddedCoshh />} />
    </Routes>
  );
};
