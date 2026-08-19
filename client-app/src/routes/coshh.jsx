import { Routes, Route } from "react-router-dom";
import { RoutedCoshhCreateModal } from "components/coshh-forms/form.jsx";

export const COSHH = () => {
  return (
    <Routes>
      <Route path="new/:reactionId" element={<RoutedCoshhCreateModal />} />

      {/* When the URL is /coshh/form/:formId,
         the Modal is gone because this route doesn't render it.
      */}
      <Route path="form/:formId/edit" element={<RoutedCoshhForm />} />
    </Routes>
  );
};
