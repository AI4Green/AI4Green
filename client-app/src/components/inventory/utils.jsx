export const EXPIRY_WARNING_DAYS = 30;

/** Whole days from today until an ISO date (negative if past), or null. */
export function daysUntil(iso) {
  if (!iso) return null;
  const today = new Date();
  today.setHours(0, 0, 0, 0);
  const target = new Date(`${iso}T00:00:00`);
  return Math.round((target - today) / 86400000);
}

export function formatDate(iso) {
  if (!iso) return "-";
  return new Date(`${iso}T00:00:00`).toLocaleDateString("en-GB", {
    day: "2-digit",
    month: "short",
    year: "numeric",
  });
}

/** 1.0000 -> "1", 0.3800 -> "0.38" */
export function formatQuantity(q) {
  return parseFloat(Number(q).toFixed(4)).toString();
}

/** Display state for a single container. */
export function itemState(item) {
  const days = daysUntil(item.expiryDate);
  if (days !== null && days < 0) return { label: "Expired", tone: "exp" };
  if (days !== null && days <= EXPIRY_WARNING_DAYS) {
    return {
      label: `Expires in ${days} day${days === 1 ? "" : "s"}`,
      tone: "exp",
    };
  }
  return item.dateOpened
    ? { label: "Opened", tone: "ok" }
    : { label: "Sealed", tone: "ok" };
}

/**
 * Roll a WorkgroupCompound's active containers up into the figures shown on
 * the compound row. Totals are summed per unit. There is no unit conversion
 * yet, so mixed units (mL and L) are shown side by side rather than merged.
 */
export function summarise(wc) {
  const active = wc.items.filter((i) => !i.isArchived);

  const totals = {};
  active.forEach((i) => {
    totals[i.unit] = (totals[i.unit] || 0) + Number(i.quantity);
  });
  const totalLabel =
    Object.entries(totals)
      .map(([unit, qty]) => `${formatQuantity(qty)} ${unit}`)
      .join(" + ") || "0";

  const locations = [...new Set(active.map((i) => i.location))];

  const expired = active.filter((i) => {
    const d = daysUntil(i.expiryDate);
    return d !== null && d < 0;
  }).length;
  const expiring = active.filter((i) => {
    const d = daysUntil(i.expiryDate);
    return d !== null && d >= 0 && d <= EXPIRY_WARNING_DAYS;
  }).length;

  const belowMin =
    wc.minStockThreshold != null &&
    (totals[wc.thresholdUnit] ?? 0) < Number(wc.minStockThreshold);

  let status;
  if (active.length === 0) status = { tone: "low", label: "Out of stock" };
  else if (expired > 0) status = { tone: "exp", label: `${expired} expired` };
  else if (expiring > 0)
    status = { tone: "exp", label: `${expiring} expiring soon` };
  else if (belowMin) {
    status = {
      tone: "low",
      label: `Low stock, minimum ${formatQuantity(wc.minStockThreshold)} ${
        wc.thresholdUnit
      }`,
    };
  } else status = { tone: "ok", label: "In stock" };

  return {
    totalLabel,
    containerCount: active.length,
    locations,
    expired,
    expiring,
    belowMin,
    status,
  };
}
