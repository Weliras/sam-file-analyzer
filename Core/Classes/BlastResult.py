from Core.Classes.SamRecord import SamRecord


class BlastResult:
    def __init__(self, sam_record: SamRecord, status_of_result: str, rid: str, no: int, result: str or None):
        self.sam_record = sam_record
        self.status_of_result = status_of_result
        self.rid = rid
        self.no = no
        self.result = result
