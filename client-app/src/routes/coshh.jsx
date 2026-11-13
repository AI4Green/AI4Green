import { Routes, Route } from "react-router-dom";
import {
  RoutedCoshhCreateModal,
  EmbeddedCoshh,
} from "components/coshh-forms/form.jsx";
import { CoshhInstanceList } from "pages/coshh-forms";

export const COSHH = () => {
  return (
    <Routes>
      <Route path="" element={<CoshhInstanceList />} />

      <Route path=":formId" element={<EmbeddedCoshh />} />
    </Routes>
  );
};
