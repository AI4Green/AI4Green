import { Routes, Route } from "react-router-dom";
import {
  RoutedCoshhCreateModal,
  EmbeddedCoshh,
} from "components/coshh-forms/form.jsx";

export const COSHH = () => {
  return (
    <Routes>
      <Route path="new" element={<RoutedCoshhCreateModal />} />

      <Route path="form/:formId/edit" element={<EmbeddedCoshh />} />
    </Routes>
  );
};
