import { Routes, Route } from "react-router-dom";

export const COSHH = () => {
  return (
    <Routes>
      <Route path="new/:reactionId" element={<CoshhCreateModal />} />

      {/* When the URL is /coshh/form/:formId,
         the Modal is gone because this route doesn't render it.
      */}
      <Route path="form/:formId/edit" element={<RoutedCoshhForm />} />
    </Routes>
  );
};
