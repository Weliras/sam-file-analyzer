import sys
import threading
import traceback
import xml.etree.ElementTree as ET

import requests
import time


from Core.Classes.BlastResult import BlastResult
from Core.Classes.SamRecord import SamRecord


blast_results = []

states = ["NO_HIT", "UNKNOWN", "FAILED", "UNEXPECTED", "SUCCESS"]

class BlastApi:
    @staticmethod
    def send_multiple_queries(program, database, sam_records, time_limit):  # type:(str, str, list[SamRecord], int) -> list[BlastResults]
        blast_results.clear()

        if len(sam_records) <= 0:
            return []

        threads = []

        start = time.time()
        for no, sam_record in enumerate(sam_records):
            t = threading.Thread(target=BlastApi.send_query, args=(program, database, sam_record, no, time_limit))
            threads.append(t)
            t.start()
        for t in threads:
            t.join()

        end = time.time()
        print(f"Time of Blast Api work: {end - start}")
        return blast_results
    @staticmethod
    def send_query(program, database, sam_record, no, time_limit):  # type:(str, str, SamRecord, int, int) -> bool
        """
        Methode which connects to blast.api and saves results from queried seq. to a file.
        :param no: number of record.
        :param program: blastn, blastp, blastx, tblastn, tblastx
        :param database: nr, cdd, nt, ...
        :param sam_record: Record that contain sequence of nucleotids ACGT
        :return: True/False if the query retrieved some hits or not
        """
        try:
            seq = sam_record.SEQ
            url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params = {"CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": seq}
            headers = {'Content-type': 'application/x-www-form-urlencoded'}
            query_request = requests.post(url=url_base, headers=headers, params=params)

            print(f"[Record no. {no}.]:{query_request.status_code}")

            # parse out the estimated time to completion
            index_rtoe = query_request.content.index(str.encode("RTOE ="))
            rtoe = int(query_request.content[index_rtoe:].decode().split("\n")[0].split("=")[1].lstrip())

            # parse out the request id
            index_rid = query_request.content.index(str.encode("RID ="))
            rid = query_request.content[index_rid:].decode().split("\n")[0].split("=")[1].lstrip()

            start = time.time()

            # wait for search to complete
            print(f"[Record no. {no}.]: Time to search {rtoe}")
            time.sleep(rtoe)

            # poll for results
            while True:
                time.sleep(60)

                url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
                params.clear()
                params = {"CMD" : "Get", "FORMAT_OBJECT" : "SearchInfo", "RID" : rid}
                result_ready_request = requests.get(url=url_base, params=params)


                index_status = result_ready_request.content.index(str.encode("Status="))
                status = result_ready_request.content[index_status:].decode().split("\n")[0].split("=")[1].lstrip()
                print(f"[Record no. {no}.]:{result_ready_request.status_code}, status: {status}")

                if status == "WAITING":
                    print(f"[Record no. {no}.]: Searching...")
                    spent_time = time.time()
                    # check if not waiting too long
                    if spent_time - start >= time_limit:
                        print(f"[Record no. {no}.]: TIME_LIMIT_EXCEEDED. Time of search: {spent_time - start}")
                        blast_results.append(BlastResult(sam_record, "TIME_LIMIT_EXCEEDED", rid, no, None))

                        url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
                        params.clear()
                        params = {"CMD": "Delete", "RID": rid}
                        result_ready_request = requests.get(url=url_base, params=params)

                        raise SystemError(
                            f"[Record no. {no}.]: Search {rid} took to long and was terminated")
                    continue
                elif status == "FAILED":
                    end = time.time()
                    print(f"[Record no. {no}.]: FAILED. Time of search: {end - start}")
                    blast_results.append(BlastResult(sam_record, "FAILED", rid, no, None))
                    raise SystemError(f"[Record no. {no}.]: Search {rid} failed; please report to blast-help\\@ncbi.nlm.nih.gov.")
                elif status == "UNKNOWN":
                    end = time.time()
                    print(f"[Record no. {no}.]: UNKNOWN. Time of search: {end - start}")
                    blast_results.append(BlastResult(sam_record, "UNKNOWN", rid, no, None))
                    raise SystemError(f"[Record no. {no}.]: Search {rid} expired.")
                elif status == "READY":
                    end = time.time()
                    print(f"[Record no. {no}.]: READY. Time of search: {end-start}")

                    index_hits = result_ready_request.content.index(str.encode("ThereAreHits="))
                    there_are_hits = result_ready_request.content[index_hits:].decode().split("\n")[0].split("=")[1].lstrip()
                    if there_are_hits == "yes":
                        print(f"[Record no. {no}.]: Search is complete, retrieving results...")
                        break
                    else:
                        print(f"[Record no. {no}.]: No hits found.")
                        blast_results.append(BlastResult(sam_record, "NO_HITS", rid, no, None))
                        return False
                else:
                    end = time.time()
                    print(f"[Record no. {no}.]: UNEXPECTED. Time of search: {end - start}")
                    blast_results.append(BlastResult(sam_record, "UNEXPECTED", rid, no, None))
                    raise SystemError(f"[Record no. {no}.]: Something unexpected happened during search. Please try it later.")

            # retrieve and display results
            url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params.clear()
            params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
            result_request = requests.get(url=url_base, params=params)
            print(f"[Record no. {no}.]:{result_request.status_code}")

            blast_results.append(BlastResult(sam_record, "SUCCESS", rid, no, result_request.content.decode()))

            #with open(f'results_{rid}.txt', 'w') as save_file:
            #    save_file.write(result_request.content.decode())

            return True

        except Exception as e:
            print(f"Oops, something bad happened : {repr(e)}")
            return False

    @staticmethod
    def analyze_results_from_blast(blast_result_list):  # type:(list[BlastResult]) -> list[str]

        found_viruses_with_count = {}

        try:
            for res in blast_result_list:
                if res.status_of_result == "SUCCESS":
                    tree = ET.ElementTree(ET.fromstring(res.result))
                    #tree = ET.parse(source="results_ZSSM6429016.txt")
                    root = tree.getroot()
                    iterations = root.find("BlastOutput_iterations").findall("Iteration")
                    for iteration in iterations:
                        iteration_hits = iteration.findall("Iteration_hits")
                        for iteration_hit in iteration_hits:
                            hits = iteration_hit.findall("Hit")
                            for hit in hits:
                                hit_def = hit.find("Hit_def")
                                if hit_def.text in found_viruses_with_count.keys():
                                    found_viruses_with_count[hit_def.text] += 1
                                else:
                                    found_viruses_with_count[hit_def.text] = 1

            interesting_results = []
            for hit_def, count in found_viruses_with_count.items():
                if count > 1:
                    interesting_results.append(hit_def)
            return interesting_results
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
            return []








